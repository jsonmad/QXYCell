# QuPath inputs, annotations, and thresholds

QXYCell reads cell measurements and optional spatial assets from a QuPath
project folder. Complete the [QuPath preparation guide](qupath_preparation.md)
([PDF](QXYCell_QuPath_Preparation_Guide.pdf)) before assembling these inputs.
Keep every exported input somewhere inside the QuPath project folder and pass
that one folder to QXYCell.

## Input files

The cell measurement table is the only unconditional input. Threshold
definitions and GeoJSON files are required only for the corresponding
downstream features.

| Input | QuPath 0.7 menu path | Requirements and use |
|---|---|---|
| **Annotation GeoJSON** *(optional)* | Select the annotation objects, then choose **File > Export objects as GeoJSON** | Export a GeoJSON `FeatureCollection` with measurements excluded. Regular labels become boolean `annotation__<label>` columns. Labels containing `Sample` define sample boundaries; labels containing a removal identifier such as `Ignore` identify regions whose cells can be excluded. |
| **Cell segmentation GeoJSON** *(optional)* | Select the cell detection objects, then choose **File > Export objects as GeoJSON** | Export a GeoJSON `FeatureCollection` with measurements excluded and QuPath Object IDs preserved. This supplies cell-boundary geometry for spatial analysis. |
| **Cell measurement table** | **Measure > Export measurements** | Select the required project images, set **Export type** to cells, and save as `measurements.csv` or `measurements.tsv`. One table may contain measurements from multiple images. |
| **Object classifier JSONs** *(alternative threshold source)* | **Classify > Object classification > Create single measurement classifier** | Create and save one simple classifier per marker. QuPath stores these under `classifiers/object_classifiers/*.json`. |
| **Threshold TSV/CSV** *(required for table-based thresholding)* | Not created in QuPath | Generate a template with `qxy.generate_threshold_table(...)`, then review and fill the per-image thresholds. A table is not required when applying object-classifier JSON thresholds directly. |

## Measurement requirements

Required columns are:

- `Image`
- `Object ID`
- `Centroid X µm`
- `Centroid Y µm`

The known encoding variants `Centroid X ¬µm` and `Centroid Y ¬µm` are accepted
and normalized automatically. Centroid columns are stored in `adata.obs` as
`Xµm` and `Yµm`.

An optional QuPath measurement column named exactly `TMA Core` is preserved and
converted into categorical `adata.obs["CoreID"]`. QXYCell does not infer
`CoreID` from any other source.

## Pixel calibration

QuPath centroid measurements are already in micrometres. GeoJSON geometry is
stored in full-resolution pixel coordinates and is scaled by
`qxy.add_annotations(..., pixel_size_um=...)`.

```python
qxy.add_annotations(adata, pixel_size_um=0.28)
```

The default is 0.28 µm. Use the verified image pixel size when it differs:

```python
qxy.add_annotations(adata, pixel_size_um=0.325)
```

The value must be positive and finite. QXYCell supports square pixels only; do
not average unequal pixel width and height values.

## Annotation behavior

`qxy.add_annotations()` imports or refreshes annotation membership and optional
cell polygons:

```python
qxy.add_annotations(adata, pixel_size_um=0.28)
```

The check report lists every annotation label, its feature count, and the
AnnData column it will populate. After import, the source-to-destination audit
is stored in `adata.uns["qxycell"]["annotation_assignments"]`, written to
`run.log`, and exported as `tables/annotation_assignments.csv`.

### Define samples

Include `Sample` in the name of each annotation that defines a sample boundary.
Matching is case-insensitive. QXYCell combines these annotations into one
categorical `adata.obs["Sample"]` column:

```python
qxy.add_annotations(adata, pixel_size_um=0.28)
adata.obs["Sample"].value_counts(dropna=False)
```

Cells inside exactly one sample annotation receive the complete annotation
name. Cells inside more than one sample annotation are labelled `Ambiguous`
and a warning is emitted. The assignment summary is stored in
`adata.uns["qxycell_sample_annotations"]`.

### Remove cells from artifact regions

Use a shared identifier such as `Ignore` in annotations drawn around tissue
artifacts, folds, debris, edge artifacts, or staining artifacts. Pass the same
identifier to `remove_cells`:

```python
adata = qxy.remove_cells(adata, remove_cells="ignore")
```

Matching is case-insensitive, so `ignore` matches names such as
`Ignore_fold`, `Ignore_edge`, and `Ignore_staining`. Because removal deletes
rows from AnnData, changed ignore polygons require reimporting measurements,
refreshing annotations, and then running removal again.

## Threshold sources

Thresholding is explicit. Choose one source for each run:

```python
# Use only QuPath single-object-classifier JSONs.
qxy.threshold_from_classifiers(adata)

# Or use only one reviewed threshold table.
qxy.threshold_from_table(adata, "/path/to/thresholds.tsv")
```

Classifier-only thresholding ignores threshold tables. Table-only thresholding
uses the named file and does not fall back to classifier JSON. Each successful
classifier run saves or replaces the applied values in
`thresholds/classifier_thresholds.tsv`.

`qxy.import_cells()` and `qxy.add_annotations()` do not apply thresholds.
`qxy.check()` reports the active threshold-definition source but does not apply
thresholds, run cell typing, or generate an LLM prompt.

### Generate and review a threshold table

```python
threshold_path = qxy.generate_threshold_table(
    "/path/to/qupath_project",
    output_dir="/path/to/outputs",
)
```

Generated tables are written under the active output folder in `thresholds/`.
Rows correspond to measurement columns containing `mean` or `median`; each
image has its own threshold column.

```text
marker    measurement_column    sample_A.tif    sample_B.tif
CD3       CD3: Mean             0.42            0.38
CD8       CD8: Median           0.31            0.29
```

The wide-table parser matches `adata.obs["Image"]` exactly. Fill one threshold
for every image column, then pass that exact file to
`qxy.threshold_from_table()`. If thresholds change, rerun thresholding and then
rerun `qxy.celltype()` to recreate threshold-dependent results.

### Refine Stage 3A thresholds with Stage 3B

A successful Stage 3A classifier run writes the thresholds it applied to
`thresholds/classifier_thresholds.tsv` in the active output folder. This table
can be used as the starting point for Stage 3B. Copy or rename it to a manual
filename such as `manual_thresholds.tsv` before editing it, because another
Stage 3A run will replace `classifier_thresholds.tsv`.

After reviewing and changing the copied values, apply that exact table with
Stage 3B:

```python
manual_threshold_table = (
    "/path/to/output/thresholds/manual_thresholds.tsv"
)
qxy.threshold_from_table(adata, manual_threshold_table)
```

Stage 3B ignores the object-classifier JSON files and uses only the named
table. Rerunning thresholding replaces the active marker-positivity columns and
marks the prompt, cell types, and post-analysis outputs as stale, so rerun
Stages 4 onward.

Classifier-derived and generated threshold tables contain one threshold column
for every unique `Image` value discovered across all measurement tables.
Image-specific classifier values populate their matching columns; a global
classifier value is repeated across all image columns. Column headings must
match `adata.obs["Image"]` exactly. Review every row and enter a numeric value
for every image being thresholded before running Stage 3B.

Recognized manual names include `thresholds.tsv`, `thresholds.csv`,
`manual_thresholds.*`, `marker_thresholds.*`, `qxycell_thresholds.*`,
`classifier_thresholds.*`, and timestamped `thresholds_*.tsv` or
`thresholds_*.csv`. When several are discovered, QXYCell prefers the most
recent timestamped file and reports which source is active.

### Conflicting classifiers

When multiple classifier JSONs define different thresholds for the same
measurement and image scope, QXYCell does not choose one by filename order.
`qxy.check()` reports every candidate and writes
`tables/classifier_conflicts.csv`. Generated tables mark the row with
`classifier_conflict=True`, preserve the candidate sources and values, and
leave image thresholds blank for review.

Measurement and annotation import can proceed, but classifier-only
thresholding refuses conflicting definitions and table-only thresholding
refuses unresolved conflict rows.

### Compartment-aware marker names

Classifier marker names may include `-nuc`, `-mem`, and `-cyto`, while cell
measurements retain the unsuffixed marker name. For example, a nucleus
classifier named `CD3` produces `CD3-nuc`. AnnData objects imported before
compartment-aware mapping was introduced must be reimported before
thresholding.

## TMA core identity

TMA core identity enters QXYCell through the cell measurement export:

1. Create and label the TMA grid in QuPath.
2. Confirm detected cells are associated with the intended cores.
3. Choose **Measure > Export measurements** and set **Export type** to cells.
4. Confirm the export includes the exact column name `TMA Core`.
5. Run `qxy.import_cells()`.

```python
adata = qxy.import_cells(project_dir)
adata.obs[["TMA Core", "CoreID"]].head()
adata.obs["CoreID"].value_counts()
```

The assignment summary is stored in
`adata.uns["qxycell_core_ids_from_measurements"]` and
`adata.uns["qxycell"]["measurement_core_assignment"]`. Without a `TMA Core`
column, QXYCell reports zero CoreIDs and does not add `CoreID`.

## Next steps

- Follow the [staged QXYCell workflow](running_qxycell.md).
- Add experimental fields with the [sample-metadata guide](metadata.md).
- Create reviewed assignments with the [cell-typing guide](cell_typing.md).
- Analyse local composition with the
  [cellular-neighbourhood guide](cellular_neighbourhoods.md).
- See the complete [function reference](QXYCell_function_reference.html).
