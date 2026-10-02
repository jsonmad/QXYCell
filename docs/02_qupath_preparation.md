# Prepare QuPath data for QXYCell

**Applies to:** QuPath 0.7.0 and QXYCell 0.1

This guide covers the QuPath actions and exports for QXYCell.

## What QXYCell needs

Keep exported inputs inside one QuPath project folder and pass that folder to QXYCell. It searches subfolders recursively. Keep QXYCell outputs beside—not inside—the project folder.

| Asset | When needed | QXYCell requirement |
|---|---|---|
| Cell measurements (`.csv` or `.tsv`) | Every run | Filename contains `measurement`, or is `detections.csv` / `detections.tsv` |
| Annotation GeoJSON | Samples, regions, or exclusions | One file per image; filename stem matches the measurement `Image` value |
| Cell GeoJSON | Cell boundaries | QuPath cell objects with measurement Object IDs; name `<image-stem>-cells.geojson` |
| Single-measurement classifier JSON | Classifier-derived thresholds | Simple classifier JSON saved below the project folder |
| Reviewed threshold table | Table-derived thresholds | Reviewed per-image TSV/CSV |

The only unconditional input is the cell measurement table. The others enable the associated QXYCell features.

## Standard route

### 1. Verify image

Open each project image and confirm the intended fluorescence series, expected
channels, and no obvious channel misregistration.

### 2. Create annotations

Use a QuPath drawing tool to create the analysis regions, give each exported
annotation a meaningful classification or name, then save the image data.

- Include `sample` anywhere in a sample-boundary name, for example
  `Sample_01` or `sample_tumour`. QXYCell imports the complete annotation name
  into a single reference column as `adata.obs["Sample"]`.
- Use unique, non-overlapping sample annotations. Cells in more than one sample
  annotation become `Ambiguous` and are reported as conflicts.
- Use a shared text string for regions to exclude, for example `Ignore_fold`
  and `Ignore_edge`. Later use the same word with
  `qxy.remove_cells(adata, remove_cells="ignore")`.
- Other annotations become boolean `annotation__<safe_label>` columns.

**Important:** For a TMA, create and label the grid before cell detection. QXYCell derives
`CoreID` only from QuPath's measurement-table column named exactly `TMA Core`,
not from annotation GeoJSON.

### 3. Segment and review cells

Run your established segmentation workflow on the required analysis regions.
InstanSeg is supported but not required; QXYCell needs QuPath cell objects,
measurements, and stable Object IDs.

For InstanSeg, open **Extensions > InstanSeg**, test an appropriate fluorescence
model and input channels on a representative region, enable **Make
measurements**, and inspect boundaries across dim, bright, crowded, sparse,
edge, and artefact regions before running the accepted settings at full scale.

Open **Measure > Show detection measurements** and confirm cells include:

- `Object ID`
- `Centroid X µm`
- `Centroid Y µm`
- intended cell, cytoplasm, or nucleus marker measurements
- `TMA Core` when TMA identity is required

QXYCell imports marker columns whose names contain `mean` or `median`
(case-insensitive). A successful segmentation run is not evidence that the
boundaries are biologically appropriate; review remains required.

### 4. Export measurements

1. Save the QuPath image data.
2. Choose **Measure > Export measurements**.
3. Select the required project images and set **Export type** to cells.
4. Export all columns unless a reduced selection has been checked against the
   required fields and marker measurements.
5. Save `measurements.csv` or `measurements.tsv` below the project folder,
   for example in `qxycell_input/`.

The required columns are exactly:

```text
Image
Object ID
Centroid X µm
Centroid Y µm
```

For a TMA, retain the exact `TMA Core` column. One export may contain several
images. Do not resave large tables in spreadsheet software: it can alter
headers, identifiers, or row counts.

### 5. Export annotations and import

For each annotated image:

1. Select annotation objects only.
2. Choose **File > Object data… > Export as GeoJSON** (or, in some QuPath 0.7
   installations, **File > Export objects as GeoJSON**).
3. Export a GeoJSON `FeatureCollection` without measurements and save it under
   the project folder.
4. Name the file using the image stem, removing only `.ome` and the image
   extension.

| Measurement `Image` value | Annotation filename |
|---|---|
| `slide01.ome.tif` | `slide01.geojson` |
| `region_A.ome.tiff` | `region_A.geojson` |
| `sample-3.tif` | `sample-3.geojson` |

`slide01-annotations.geojson` does not match `slide01.ome.tif`. QXYCell uses
this stem match to assign annotation geometry to the correct cell rows.

Then import cells and GeoJSON:

```python
adata = qxy.import_cells(project_dir)
qxy.add_annotations(adata)
```

Inspect a spatial overlay before downstream analysis.

## Calibration, annotations, and thresholds

### Pixel calibration

In QuPath's **Image** tab, check that pixel width and height are equal. QuPath
centroids are already in µm, but GeoJSON uses full-resolution pixels. QXYCell
uses `0.28` µm/pixel by default; supply the verified value only when it differs:

```python
qxy.add_annotations(adata, pixel_size_um=0.325)
```

Do not average unequal pixel dimensions.

### Annotation behavior

Annotations containing `Sample` create `adata.obs["Sample"]`; other labels
create `annotation__<label>` columns. Use `Ignore` in artifact-region labels
and remove those cells with:

```python
adata = qxy.remove_cells(adata, remove_cells="ignore")
```

Changed Ignore polygons require reimporting measurements, refreshing
annotations, and removing cells again.

### Threshold tables and refinement

Apply classifier JSON thresholds directly with
`qxy.threshold_from_classifiers(adata)`, or apply a reviewed table with
`qxy.threshold_from_table(adata, threshold_file)`. A classifier run writes
`thresholds/classifier_thresholds.tsv` with one threshold column per image.
Copy or rename that table before manually adjusting values, then apply the
reviewed copy with `qxy.threshold_from_table()`.

To begin from a fresh table instead, run:

```python
threshold_file = qxy.generate_threshold_table(project_dir)
# Review every marker and image value, then:
qxy.threshold_from_table(adata, threshold_file)
```

Classifier-only thresholding refuses conflicting JSON definitions. Generated
tables preserve conflict information and require a numeric value for every
affected image before table thresholding can run.

## Optional assets

### TMA core identity

Create and label the QuPath TMA grid before detection. In **Measure > Export
measurements**, export cells and retain the exact `TMA Core` field.
`qxy.import_cells()` preserves it and creates categorical `CoreID`; without it,
QXYCell does not infer core identity. See QuPath's
[TMA grid commands](https://qupath.readthedocs.io/en/stable/docs/reference/commands.html#tma).

### Threshold classifiers

After segmentation, create one reviewed starting threshold per marker:

1. Choose **Classify > Object classification > Create single measurement classifier**.
2. Filter to cells and select an appropriate mean or median measurement.
3. Set the below/above-threshold classes, enable live preview, and review
   representative images and tissue conditions.
4. Save a short, unique marker-based classifier name below the project folder.

QXYCell reads simple single-measurement classifier JSONs. Composite or malformed
classifiers are reported but not converted into threshold rows. Classifier
thresholds can be applied directly, or copied from the generated
`classifier_thresholds.tsv` table, manually refined, and applied as a reviewed
table.

### Cell-boundary GeoJSON

Cell GeoJSON is optional, but required for `cell_polygon_wkt` and cell-boundary
plots. For each image:

1. Select cell detection objects, not annotations.
2. Choose the GeoJSON export command and export a `FeatureCollection` without
   measurements.
3. Preserve QuPath Object IDs and save as `<image-stem>-cells.geojson`, for
   example `qxycell_input/cells/slide01-cells.geojson`.

QXYCell distinguishes the two exports through QuPath's GeoJSON `objectType`:
`annotation` features are assigned by their matching image-stem filename,
while `cell` features are matched to measurements by Object ID. Do not
regenerate Object IDs between the measurement and cell-GeoJSON exports. For
large images, export one image at a time and sanity-check the feature count.

### Suggested input layout

```text
qupath_project/
|-- project.qpproj
|-- data/
|-- classifiers/
|   `-- object_classifiers/
|       |-- CD3.json
|       `-- PanCK.json
`-- qxycell_input/
    |-- measurements.tsv
    |-- annotations/
    |   `-- slide01.geojson
    `-- cells/
        `-- slide01-cells.geojson
```

Source OME-TIFFs can remain outside the project when project links are valid;
the exported inputs above must remain inside it.

## Troubleshooting

| Symptom | Recovery |
|---|---|
| No measurement file is found | Rename it to include `measurement`, or use `detections.csv` / `detections.tsv`. |
| Required columns are missing | Re-export cells with `Image`, `Object ID`, `Centroid X µm`, and `Centroid Y µm`; do not edit headers. |
| Annotations are found but not assigned | Match each GeoJSON filename stem to `Image`. |
| Cell polygons are missing | Export cell objects rather than annotations and preserve the matching Object IDs. |

## Final checklist

- [ ] Intended image series and channels verified.
- [ ] Segmentation reviewed across representative regions.
- [ ] Measurement export contains every required column and, for TMA, `TMA Core`.
- [ ] Annotation filenames match their image stems.
- [ ] Cell GeoJSON is named `<image-stem>-cells.geojson` and retains the
  measurement-export Object IDs, when used.
- [ ] Classifier JSONs are simple single-measurement classifiers, when used.
- [ ] Spatial overlay alignment reviewed after import.

## Official references

- [QuPath 0.7.0 releases](https://github.com/qupath/qupath/releases)
- [QuPath projects](https://qupath.readthedocs.io/en/stable/docs/tutorials/projects.html)
- [QuPath image concepts and pixel calibration](https://qupath.readthedocs.io/en/stable/docs/concepts/images.html)
- [QuPath InstanSeg workflow](https://qupath.readthedocs.io/en/stable/docs/deep/instanseg.html)
- [QuPath measurement export](https://qupath.readthedocs.io/en/stable/docs/tutorials/exporting_measurements.html)
- [QuPath GeoJSON export](https://qupath.readthedocs.io/en/stable/docs/advanced/exporting_annotations.html)
- [QuPath 0.7.0 command reference](https://qupath.readthedocs.io/en/stable/docs/reference/commands.html)
