# AnnData structure and outputs

QXYCell stores measurements, spatial coordinates, annotations, analysis
results, palettes, output paths, and provenance in one AnnData object. Every
successful staged operation refreshes the active H5AD checkpoint and relevant
table exports.

## AnnData layout

| Location | Added by | Contents |
|---|---|---|
| `adata.obs["Image"]` | `qxy.import_cells()` | QuPath image name per cell |
| `adata.obs["Xµm"]`, `adata.obs["Yµm"]` | `qxy.import_cells()` | Cell centroid coordinates in micrometres |
| `adata.obs["<marker>_pos"]` | threshold stages | Boolean marker positivity columns |
| `adata.obs["annotation__<label>"]` | `qxy.add_annotations()` | Boolean annotation membership columns |
| `adata.obs["cell_polygon_wkt"]` | `qxy.add_annotations()` / `qxy.load_cell_polygons()` | Optional cell segmentation geometry |
| `adata.obs["Sample"]` | `qxy.add_annotations()` / `qxy.assign_samples()` | Sample label from annotations containing `Sample` |
| `adata.obs["TMA Core"]` | `qxy.import_cells()` | QuPath measurement-table TMA core label |
| `adata.obs["CoreID"]` | `qxy.import_cells()` / `qxy.assign_core_ids_from_measurements()` | Categorical CoreID derived only from `TMA Core` |
| `adata.obs["celltype"]` | `qxy.celltype()` | Assigned cell type |
| `adata.obs["cn"]` | `qxy.cn_kmeans()` / `qxy.cn_name()` | Neighbourhood cluster label and optional descriptive name |
| `adata.obs[*metadata columns*]` | `qxy.add_metadata()` | Sample metadata broadcast to cells |
| `adata.X` | `qxy.import_cells()` | Marker intensity matrix, cells × markers |
| `adata.var` | `qxy.import_cells()` | Marker names and metadata |
| `adata.obsm["spatial"]` | `qxy.import_cells()` | Cell centroid x/y coordinates in micrometres |
| `adata.obsm["cn_profile"]` | `qxy.cn_knn()` | Per-cell local cell-type composition |

## Provenance and summaries

| Location | Contents |
|---|---|
| `adata.uns["qxycell"]` | Project and output paths, timestamps, stage status, provenance, palettes, and active H5AD path |
| `adata.uns["qxycell_annotation_labels"]` | Annotation class to `annotation__<label>` column mapping |
| `adata.uns["qxycell_thresholding"]` | Threshold source and positivity-column summary |
| `adata.uns["qxycell_sample_annotations"]` | Sample assignment summary |
| `adata.uns["qxycell_core_ids_from_measurements"]` | Measurement-derived CoreID summary |
| `adata.uns["qxycell_dataset_summary"]` | Paths to dataset-summary tables and HTML report |
| `adata.uns["qxycell_sample_metadata"]` | Metadata matching summary |
| `adata.uns["qxycell_celltyping"]` | Cell-typing rule summary |
| `adata.uns["cn"]` | Neighbourhood parameters, cell-type list, and label map |

`run.log` records the same staged provenance, including threshold and YAML
sources and whether thresholds, cell typing, or prompt generation were
performed.

## Active output folder

`qxy.import_cells()` creates a timestamped sibling run folder and stores its
path in `adata.uns["qxycell"]["output_dir"]`:

```text
qupath_project_run_YYMMDD_HHMM/
  h5ad/
    qxycell.h5ad
  tables/
    cells_obs.csv
    markers_var.csv
    annotation_assignments.csv
  thresholds/
    classifier_thresholds.tsv
  celltype/
    current_prompt.txt
  dataset_summary/
    dataset_summary.html
    dataset_overview.tsv
  cn/
    cn_labels.csv
  plots/
  run.log
```

Only outputs created by the stages that were run will be present. Downstream
functions reuse the active folder when `output_dir` is omitted.

## Dataset summary

Generate descriptive tables and an HTML report for the cells, markers,
annotations, samples, and cell types already present in AnnData:

```python
summary = qxy.dataset_summary(adata, sample_col="Sample")
```

Use `sample_col="Image"` when sample annotations are unavailable. The report is
written under `dataset_summary/` in the active output folder. Depending on the
fields available in AnnData, outputs include:

- `dataset_summary.html`;
- `dataset_overview.tsv`;
- `cells_per_sample.tsv`;
- `celltype_counts.tsv` and `celltypes_per_sample.tsv`;
- `marker_positivity.tsv` and `marker_positivity_by_sample.tsv`;
- `annotation_counts.tsv`; and
- `ignore_cells.tsv`.

The returned dictionary contains the generated tables and paths. Paths are
also stored in `adata.uns["qxycell_dataset_summary"]`. Save the AnnData object
afterward if that provenance should be retained in the H5AD:

```python
qxy.save(adata)
```

This is descriptive reporting. It does not validate image quality,
segmentation, staining, thresholds, batch effects, cell-type accuracy, or
spatial alignment.

## Save

`qxy.save(adata)` writes a compressed H5AD containing the full analysis state:

```python
# Update the active run checkpoint.
qxy.save(adata)

# Save under a chosen output folder.
qxy.save(adata, output_dir="my_qxy_output")

# Save to an exact file path.
qxy.save(adata, path="my_analysis.h5ad")
```

The default path is:

```text
qupath_project_run_YYMMDD_HHMM/h5ad/qxycell.h5ad
```

If `adata.uns["qxycell"]["h5ad_path"]` already exists, `qxy.save(adata)`
updates that file.

## Load

```python
# Load the most recent QXYCell save discoverable from the current location.
adata = qxy.load_latest()

# Load one exact H5AD file.
adata = qxy.load("path/to/qxycell.h5ad")
```

Reloading restores marker values, cell metadata, spatial coordinates,
annotation summaries, palettes, neighbourhood results, and staged provenance
already saved in the object.

## CSV and TSV review files

The active H5AD is the analysis checkpoint. CSV and TSV outputs are review and
interchange files that mirror important parts of the object:

- `tables/cells_obs.csv` mirrors per-cell `adata.obs` fields.
- `tables/markers_var.csv` mirrors marker-level `adata.var` fields.
- `tables/annotation_assignments.csv` records annotation source, destination,
  and assigned-cell counts.
- `thresholds/classifier_thresholds.tsv` records classifier-derived thresholds
  actually applied.
- Dataset-summary and plotting functions save their own matrices and summaries
  beside generated reports and figures.

After a successful stage, use the refreshed H5AD as the source for subsequent
analysis rather than reconstructing state from these flat files.

## Related documentation

- Understand checkpoint replacement in
  [Running the staged workflow](running_qxycell.md).
- Review annotation and threshold provenance in
  [QuPath inputs, annotations, and thresholds](qupath_inputs.md).
- Add experimental fields with the [sample-metadata guide](metadata.md).
- Assign reviewed labels with the [cell-typing guide](cell_typing.md).
- Analyse local composition with the
  [cellular-neighbourhood guide](cellular_neighbourhoods.md).
- Create figures with the [plotting guide](plotting.md).
