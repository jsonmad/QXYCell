# QXYCell

<p align="center">
  <img src="assets/qxycell-icon.png" alt="QXYCell icon" width="180">
</p>

QXYCell converts single-cell measurements and spatial assets from a QuPath
project folder into analysis-ready AnnData `.h5ad` objects for cell typing,
spatial analysis, and plotting.

The resulting AnnData object can be used with downstream tools such as
[Scanpy](https://scanpy.readthedocs.io/en/stable/),
[Squidpy](https://squidpy.readthedocs.io/en/stable/), and
[scimap](https://scimap.xyz/), or incorporated into a
[SpatialData](https://spatialdata.scverse.org/en/stable/) workflow when image,
label, or shape elements are also needed.

## Contents

- [Installation](#installation)
- [Prepare data in QuPath](#prepare-data-in-qupath)
- [Quick start](#quick-start)
- [Documentation and support](#documentation-and-support)
- [QuPath inputs](#qupath-inputs)
- [Annotations](#annotations)
- [QC](#qc)
- [Metadata](#metadata)
- [Cell typing](#cell-typing)
- [Cellular neighbourhoods](#cellular-neighbourhoods)
- [Plots and heatmaps](#spatial-plots)
- [Save and load](#save-and-load)
- [Optional workflow shortcut](#optional-workflow-shortcut)

## Installation

QXYCell requires Python 3.10 or newer.

### Create the QXYCell environment

```bash
git clone https://github.com/jsonmad/QXYCell.git
cd QXYCell
conda env create -f environment.yml
conda activate qxycell
```

The environment file installs QXYCell and its scientific Python dependencies
into the new environment.

Verify the installation:

```bash
python -c "import qxycell; print('QXYCell import OK')"
qxycell --help
```

### Updating QXYCell

When a new QXYCell version is announced, update it manually from inside the
cloned repository:

```bash
conda activate qxycell
git pull
conda env update -f environment.yml --prune
```

![QXYCell workflow from multiplex tissue imaging through QuPath, QXYCell,
AnnData, spatial plots, and downstream analysis](docs/assets/qxycell_workflow.png)

> QXYCell is independent and is not affiliated with or endorsed by
> [QuPath](https://qupath.github.io/) ([GPLv3](https://github.com/qupath/qupath/blob/main/LICENSE)).

## Prepare data in QuPath

Before running QXYCell, follow the
[QuPath preparation guide](docs/qupath_preparation.md) ([PDF version](docs/QXYCell_QuPath_Preparation_Guide.pdf)).
It covers image and pixel-size verification, annotations, InstanSeg cell
segmentation, measurement and classifier export, GeoJSON export, filenames,
and the required preflight check. The workflow targets QuPath 0.7.0 and is
applicable to multiplex immunofluorescence data from any acquisition platform.

A multichannel OME-TIFF does not need to originate from COMET, but it must open
correctly in QuPath and have verified channels, dimensions, registration, and
physical pixel calibration. QXYCell supports square pixels only.

## Quick start

```python
import qxycell as qxy

# Optional: validate your QuPath project folder before running
report = qxy.check("/path/to/qupath_project")

# Stage 1: import the measurement table into AnnData
adata = qxy.import_cells("/path/to/qupath_project")

# Stage 2: add or refresh GeoJSON annotations and cell polygons
qxy.add_annotations(adata, pixel_size_um=0.28)

# Optional Stage 2b: remove cells in tissue or staining artifact regions
qxy.remove_cells(adata, remove_cells="ignore")

# Stage 3A: use QuPath classifier JSON thresholds only
qxy.threshold_from_classifiers(adata)
# Saves the applied values to thresholds/classifier_thresholds.tsv

# Alternatively, stage 3B uses one named threshold table only:
# table = qxy.generate_threshold_table("/path/to/qupath_project")
# qxy.threshold_from_table(adata, table)

# Stage 4: create the prompt used to draft celltype_logic.yaml
qxy.celltype_prompt(adata, context="Describe the tissue and expected populations")

# Pause for expert review, then save the returned YAML
# Stage 5: assign cell types from the reviewed YAML
qxy.celltype(adata, "/path/to/celltype_logic.yaml")

# Optional Stage 6: plot assigned cell types without opening GUI windows
qxy.plot_spatial(adata, category_col="celltype", show=False)
```

### Why the workflow is staged

Stage 1 creates the base measurement checkpoint. Stages 2–5 update that same
H5AD and refresh `tables/cells_obs.csv` and `tables/markers_var.csv`, so the
analysis can be reviewed and revised without reimporting the measurements. The
optional Stage 2b removes cells inside user-labelled regions containing tissue
artifacts, folds, debris, edge artifacts, or staining artifacts, then refreshes
the filtered H5AD and `tables/cells_obs.csv`.

| When an input changes | Rerun | What QXYCell replaces |
|---|---|---|
| Annotation or cell GeoJSON | `qxy.add_annotations(adata)` | Annotation, sample, and cell-polygon columns; downstream stages become stale |
| Ignore annotation polygons | Stages 1 and 2, then `qxy.remove_cells(adata, remove_cells="ignore")` | Rebuilds all cells before removing those inside tissue or staining artifact regions |
| QuPath classifier JSON thresholds | `qxy.threshold_from_classifiers(adata)` | Marker `_pos` columns and `thresholds/classifier_thresholds.tsv`; prompt, cell types, and post-analysis become stale |
| A reviewed threshold table | `qxy.threshold_from_table(adata, table)` | Marker `_pos` columns; prompt, cell types, and post-analysis become stale |
| Biological context for the prompt | `qxy.celltype_prompt(adata, context=...)` | `celltype/current_prompt.txt`; an expert-edited YAML is preserved |
| Cell-type YAML | `qxy.celltype(adata)` | `celltype`, feature, derived-feature, count, and rule-summary outputs |

The two threshold functions are deliberately separate. Classifier-only
thresholding ignores tables; table-only thresholding uses the named table and
does not fall back to classifier JSON.

Each successful classifier-only threshold run replaces
`thresholds/classifier_thresholds.tsv` in the active output folder. This table
records the classifier-derived values that were applied and can be reviewed or
passed later to `qxy.threshold_from_table()`.

For runnable files with one Python script per stage, see the
[staged workflow examples](staged_workflow/README.md). The
[interactive notebook](staged_workflow/QXYCell_staged_workflow.ipynb)
and numbered scripts implement the same staged workflow; choose the notebook
for a cell-by-cell experience or the scripts for terminal use.

When the verified QuPath pixel size differs from 0.28 µm, supply the single
square-pixel value during import:

```python
qxy.add_annotations(adata, pixel_size_um=0.325)
```

The value must be positive and finite. QuPath centroid measurements are already
in micrometres; `pixel_size_um` scales annotation and cell GeoJSON coordinates
from full-resolution pixels into micrometres. Do not average unequal pixel width
and height values—QXYCell does not support non-square pixels.

When called, `qxy.check()` writes a timestamped sibling check folder.
`qxy.import_cells()` creates the timestamped sibling run folder used by
all later stages:

```text
qupath_project_check_YYMMDD_HHMM/
qupath_project_run_YYMMDD_HHMM/
```

Downstream functions such as annotation import, the explicit threshold stages,
`qxy.add_metadata()`, `qxy.celltype()`, plotting, and `qxy.save()` reuse the
active run folder stored in `adata.uns["qxycell"]["output_dir"]` when
`output_dir` is omitted.

**The measurement and annotation stages populate AnnData as follows:**

| Location | Contents |
|---|---|
| `adata.obs` | Per-cell metadata: `Image`, `Object ID`, `Xµm`, `Yµm`, optional measurement-table `TMA Core` and derived `CoreID`, `annotation__<label>` boolean columns, `Sample` when sample annotations exist, and `cell_polygon_wkt` when cell GeoJSON is available |
| `adata.X` | Marker intensity matrix (cells × markers) |
| `adata.var` | Marker names and metadata |
| `adata.obsm["spatial"]` | Cell centroid coordinates in microns, shape `(n_cells, 2)` |
| `adata.uns["qxycell"]` | Run metadata: output directory, timestamps, file paths |
| `adata.uns["qxycell_annotation_labels"]` | Map of annotation class names to `annotation__<label>` column names |

## Documentation and support

The repository includes a [QuPath preparation guide](docs/qupath_preparation.md),
an [overview](docs/QXYCell_overview.html), a
[function reference](docs/QXYCell_function_reference.html), and
[synthetic function examples](docs/qxy_function_examples.html). The rendered
documentation site will be enabled after the repository becomes public.

Report reproducible bugs and feature requests through
[GitHub Issues](https://github.com/jsonmad/QXYCell/issues). Do not attach
identifiable patient or research data, credentials, or other sensitive material
to a public issue. Report suspected vulnerabilities privately as described in
[SECURITY.md](SECURITY.md).

QXYCell is released under the [MIT License](LICENSE).

## QuPath inputs

QXYCell reads exported assets from a QuPath project folder. Complete the
[QuPath preparation and verification workflow](docs/qupath_preparation.md)
before assembling these inputs. Keep every exported input somewhere inside the
QuPath project folder and pass that one folder to QXYCell; do not create a
second input directory. The measurement table is the
only unconditional input. Threshold definitions and GeoJSON files are needed
only for the corresponding downstream features:

- **Cell measurement table** — `measurements.csv` or `measurements.tsv` exported from QuPath. One table may contain cells from multiple images.
- **Threshold TSV/CSV** *(required for thresholding)* — the source of truth for marker positivity thresholds. Use a filled table such as `thresholds.tsv` or `thresholds_YYMMDD-HHMM.tsv`. Generated tables are written under the QXYCell output folder in `thresholds/`.
- **Object classifier JSONs** *(alternative threshold-template source)* — single-measurement classifiers saved under `classifiers/object_classifiers/*.json`. When no threshold table exists, QXYCell can convert these JSONs into a fresh timestamped table for review. Existing threshold tables remain the active source.
- **Annotation GeoJSON** *(optional)* — exported QuPath annotation polygons, with measurements excluded. Regular annotation classification/name labels become boolean `annotation__<label>` columns in `adata.obs`. Annotations with `Sample` in the label define sample boundaries and are collapsed into one `adata.obs["Sample"]` column; annotations labelled `Ignore` mark regions to exclude. Annotation labels never create or replace `CoreID` values.
- **Cell segmentation GeoJSON** *(optional)* — exported cell objects for all cells, measurements excluded. Provides geometry for spatial analysis.

Required measurement columns: `Image`, `Object ID`, `Centroid X µm`, `Centroid Y µm`.
The known encoding variants `Centroid X ¬µm` and `Centroid Y ¬µm` are accepted and
normalized automatically.
These QuPath centroid columns are stored in `adata.obs` as `Xµm` and `Yµm`.
QuPath GeoJSON geometry is stored in full-resolution pixel coordinates and is
scaled by `qxy.add_annotations(..., pixel_size_um=...)`; the default is `0.28`. Verify the
image's square-pixel calibration in QuPath before running QXYCell.
An optional QuPath measurement column named exactly `TMA Core` is preserved and
converted into categorical `adata.obs["CoreID"]` by `qxy.import_cells()`. QXYCell does
not infer `CoreID` from any other source. If the measurement
table has no `TMA Core` column, the resulting AnnData has no `CoreID` column
and the check report lists zero CoreIDs.

Threshold tables are used explicitly by `qxy.threshold_from_table()`.
Generated threshold tables are written to a `thresholds/` folder inside the
QXYCell output folder:

```text
outputs/
  h5ad/
  tables/
  thresholds/
    thresholds_YYMMDD-HHMM.tsv
  run.log
```

The table includes only measurement columns whose names contain `mean` or
`median`. Each row is one feature, and each image gets its own threshold
column. When a usable object classifier JSON maps to that measurement column,
the `marker` value is taken from the classifier JSON filename, not from the
measurement column text. For example, `aSMA.json` mapped to
`Cell: #945;SMA - TRITC: Mean` writes marker `aSMA`.

`qxy.check()` does not generate threshold tables. To generate a fresh threshold
table from object classifier JSONs explicitly, call:

```python
threshold_path = qxy.generate_threshold_table("/path/to/qupath_project")
```

The check report lists every annotation name found in GeoJSON, its feature
count, and the AnnData column it will populate. Sample-labelled annotations are
shown together as inputs to `Sample`; other labels map to
`annotation__<safe_label>` columns.

The report also separates definitions from actions. `qxy.check()` reports the
active threshold-definition source but explicitly records that thresholds,
cell typing, and LLM prompt generation were not performed. `check_report.json`
stores these flags under `processing`.

After annotation import, `qxy.add_annotations()` reports the same source-to-destination mapping with
the number of cells actually assigned. The audit is stored in
`adata.uns["qxycell"]["annotation_assignments"]`, written to `run.log`, and
exported as `tables/annotation_assignments.csv`.

`run.log` and `adata.uns["qxycell"]` also record whether thresholding and cell
typing were applied, the threshold or YAML source used, and whether an LLM
prompt was generated. Defaults are `apply_thresholds=False`,
`celltype_logic=None`, and no prompt generation. Prompt creation is always an
explicit `qxy.celltype_prompt(adata)` call.

To control where that table is written, pass the same output folder used for
the run:

```python
threshold_path = qxy.generate_threshold_table(
    "/path/to/qupath_project",
    output_dir="/path/to/outputs",
)
```

`qxy.import_cells()` and `qxy.add_annotations()` do not apply or select
thresholds. Choose `qxy.threshold_from_classifiers()` to use only QuPath
classifier JSON files, or `qxy.threshold_from_table()` to use only the named
table. Neither function silently falls back to the other source.

### Conflicting classifiers and per-image review

When multiple classifier JSONs contain different thresholds for the same
measurement and image scope, QXYCell no longer selects one by filename order.
`qxy.check()` reports every candidate and writes
`tables/classifier_conflicts.csv`. Generated threshold tables mark the row with
`classifier_conflict=True`, preserve candidate names, values, and sources, and
leave every image threshold blank for review. Measurement and annotation import
can still proceed, but classifier-only thresholding refuses conflicting
definitions and table-only thresholding refuses unresolved conflict rows.

Fill one threshold for every image column in the generated table, then reapply
thresholding using that exact file:

```python
threshold_file = output_dir / "thresholds" / "thresholds_YYMMDD-HHMM.tsv"

qxy.threshold_from_table(adata, threshold_file)
qxy.celltype(adata, "celltype_logic.yaml")
```

Different image columns may contain different values. The wide-table parser
matches `adata.obs["Image"]` exactly and applies the corresponding threshold to
each image. After editing the table again, rerun the same calls. Thresholding
removes threshold-dependent cell-type and feature columns; the subsequent
`qxy.celltype()` call recreates active
labels from the new per-image positivity calls. The active H5AD and
`tables/cells_obs.csv` are refreshed after every stage.

To manually edit thresholds, fill in or edit the image threshold columns, save
the finished file in the QuPath project folder or in the active run's
`thresholds/` subfolder, then pass that exact path to
`qxy.threshold_from_table(adata, threshold_file)`. `qxy.check()` can still be
given the same file to validate it. The compact `CheckReport` summary includes
`Threshold source: ...`, and `check_report.txt` includes
`Active threshold source: ...`.

If multiple threshold files are found, QXYCell prefers timestamped
`thresholds_*.tsv`/`.csv` files and uses the most recently modified one.
Other threshold files are ignored with a warning in the check report.

Recognized manual threshold filenames are `thresholds.tsv`, `thresholds.csv`,
`manual_thresholds.tsv`, `manual_thresholds.csv`, `marker_thresholds.tsv`,
`marker_thresholds.csv`, `qxycell_thresholds.tsv`, `qxycell_thresholds.csv`,
`classifier_thresholds.tsv`, `classifier_thresholds.csv`, and timestamped
`thresholds_*.tsv` / `thresholds_*.csv`.

```text
marker    measurement_column    sample_A.tif    sample_B.tif
CD3       CD3: Mean             0.42            0.38
CD8       CD8: Median           0.31            0.29
```

Classifier marker names include non-cellular compartments: `-nuc`, `-mem`,
and `-cyto`. Cell measurements keep the unsuffixed marker name. For example,
a nucleus classifier named `CD3` produces the classifier marker `CD3-nuc`,
while the cell measurement remains `CD3`. AnnData objects imported before this
change must be re-imported before thresholding so the threshold-to-intensity
mapping is available.

## Annotations

Remove cells inside annotations drawn around tissue artifacts, tissue folds,
debris, edge artifacts, or staining artifacts:

```python
adata = qxy.remove_cells(adata, remove_cells="ignore")
```

The match is case-insensitive. Rows inside matching annotation polygons are
removed from `adata.obs` in place.

`qxy.add_annotations()` converts annotations with `Sample` in the label into one `Sample`
column. Cells inside more than one sample annotation are labelled `Ambiguous`
and a warning is emitted. You can rerun sample assignment explicitly:

```python
sample_summary = qxy.assign_samples(adata)
adata.obs["Sample"].value_counts()
```

Adds `adata.obs["Sample"]` (string). Summary stored in `adata.uns["qxycell_sample_annotations"]`.

## QC

Generate per-sample QC tables and an HTML report:

```python
qc = qxy.qc(adata, sample_col="Image")
```

Results stored in `adata.uns["qxycell_qc"]`. HTML report and TSV tables are
written to the active output folder under `qc/`.

## Metadata

Add sample-level metadata matched on image name:

```python
qxy.add_metadata(adata, "sample_metadata.tsv")
```

Match on a custom sample column or a different key column in the metadata file:

```python
qxy.add_metadata(
    adata,
    "sample_metadata.tsv",
    sample_col="ImageID",          # column in adata.obs to match on
    metadata_sample_col="sample",  # column in the TSV (if different)
)
```

Import only selected columns:

```python
qxy.add_metadata(adata, "sample_metadata.tsv", columns=["group", "mouse_id"])
```

Each column in the metadata file is broadcast to every cell belonging to that sample and added as a new column in `adata.obs`. For example, a `group` column in the TSV becomes `adata.obs["group"]`. Summary stored in `adata.uns["qxycell_sample_metadata"]`.

## Cell typing

### How thresholds are used for cell typing

Cell typing uses marker-positive calls rather than the raw measurement values
directly. First apply thresholds using one of the two explicit threshold stages:

- `qxy.threshold_from_classifiers(adata)` reads the simple thresholds from the
  QuPath object-classifier JSON files, applies them, and saves the values used
  to `thresholds/classifier_thresholds.tsv` in the active output folder.
- `qxy.threshold_from_table(adata, threshold_file)` ignores classifier JSON
  thresholds and applies only the reviewed threshold table supplied by the
  user.

For each marker, QXYCell compares the selected QuPath measurement value with
the applicable threshold. Thresholds can be global or specific to each image.
A value greater than or equal to the threshold is recorded as marker-positive
in `adata.obs["<MARKER>_pos"]`; a value below the threshold is recorded as
marker-negative.

The cell-type YAML then refers to these positivity columns by marker name.
`positive` requires all listed markers to be positive, `negative` requires the
listed markers to be negative, and `any_positive` requires at least one listed
marker to be positive. Rules are evaluated from top to bottom, and the first
matching rule assigns `adata.obs["celltype"]`. Put specific cell types before
broader populations so broad rules do not capture them first.

After changing to compartment-aware marker names, regenerate the cell-type
prompt and update YAML marker references before applying cell types again. For
example, a rule that previously referenced `CD3` for a nucleus classifier must
reference `CD3-nuc`.

If thresholds are changed and either threshold stage is rerun, QXYCell replaces
the previous marker-positive calls and removes threshold-dependent cell-type
and feature columns. Rerun `qxy.celltype()` afterwards to recreate the cell
typing results from the updated thresholds.

Generate a first-pass cell type prompt for an LLM:

```python
prompt = qxy.celltype_prompt(adata)
print(prompt)
```

This prints the prompt, replaces `celltype/current_prompt.txt` in the active
output folder,
and returns the string. After editing and saving a YAML logic file in that
folder, apply it:

```python
summary = qxy.celltype(adata)
```

Adds `adata.obs["celltype"]` (string). Cell type assignment summary stored in `adata.uns["qxycell_celltyping"]`. QXYCell applies the most recently saved YAML in the celltype folder and prints the file path it used. Re-running either explicit threshold function removes threshold-dependent cell-type and feature columns, marks them stale in the stage record, and requires cell typing to be rerun. Re-running cell typing replaces columns owned by the previous YAML, including obsolete feature columns.

## Cellular neighbourhoods

Compute cellular neighbourhoods (CNs) using k-nearest neighbours:

```python
# Step 1: build KNN graph and compute per-cell local composition profiles
qxy.cn_knn(adata, k=15)
```

Adds `adata.obsm["cn_profile"]` — a float32 array of shape `(n_cells, n_cell_types)` where each row is the local cell type composition around that cell. Run parameters stored in `adata.uns["cn"]`.

```python
# Step 2: cluster the profiles into CN groups
qxy.cn_kmeans(adata, n_cn=10)
```

Adds `adata.obs["cn"]` (integer cluster label). Updates `adata.uns["cn"]` with clustering parameters.

Auto-assign descriptive names to CN clusters from their composition profiles:

```python
label_table = qxy.cn_name(adata)
```

Labels are derived from `adata.obsm["cn_profile"]` using a priority rule: clusters where one cell type dominates (≥ 50 %) are labelled `<type> hi`; mixed clusters show the top two contributors (e.g. `CD8 T + Macrophage`). Renames values in `adata.obs["cn"]` (string). Label map stored in `adata.uns["cn"]["label_map"]`. Summary saved to the active output folder under `cn/cn_labels.csv`.

Optionally shorten long cell type names before labelling:

```python
label_table = qxy.cn_name(adata, compaction={"CD8+PD1+LAG3+": "PD1 LAG3 CD8"})
```

Generated CN labels are Windows-safe. Path separators are rendered as ` + `,
Windows-reserved filename characters are removed or replaced, repeated `..`
segments are collapsed, and trailing spaces or periods are stripped.

## Spatial plots

Plot cell type distributions in tissue space:

```python
# Automatically use Sample when available, otherwise Image
qxy.plot_spatial(adata)

# Force one plot per QuPath image
qxy.plot_spatial(adata, sample_col="Image")

# One plot per Sample annotation, when a Sample column is available
qxy.plot_spatial(adata, sample_col="Sample")

# Combined multi-panel figure
qxy.plot_spatial(adata, combined=True, max_cols=3)

# Selected images only when explicitly using Image grouping
qxy.plot_spatial(adata, samples=["image_A.ome.tif", "image_B.ome.tif"])

# Fixed square window (microns)
qxy.plot_spatial(adata, fixed_window_um=11500)

# Cell-distribution centering instead of the default bounding-box center
qxy.plot_spatial(adata, center_method="median")

# Figure panel aspect from selected sample X/Y extent
qxy.plot_spatial(adata, auto_figsize=True)

# Use raw coordinate orientation instead of the default image-viewer y flip
qxy.plot_spatial(adata, flip_y=False)

# Plot CNs instead of cell types
qxy.plot_spatial(adata, category_col="cn")
```

When `sample_col` is omitted, `qxy.plot_spatial()` uses
`adata.obs["Sample"]` when that column contains usable labels and otherwise
falls back to `adata.obs["Image"]`. This allows datasets without QuPath Sample
annotations to plot normally. Pass `sample_col="Sample"` or
`sample_col="Image"` to force one grouping method. The `samples=` argument
always selects values from the resolved `sample_col`.

Cells with a missing value in `adata.obs[sample_col]` are excluded by default.
Use `include_missing_samples=True` only when a separate `"nan"` panel is wanted.

Choose spatial plot file formats:

```python
qxy.plot_spatial(adata, save_png=True, save_pdf=False)   # PNG only
qxy.plot_spatial(adata, save_png=False, save_pdf=True)   # PDF only
qxy.plot_spatial(adata, save_png=True, save_pdf=True)    # both
```

Use a full AnnData object as the grey underlay while plotting categories from
a filtered object:

```python
qxy.plot_spatial(
    adata_cn,
    underlay_adata=adata,
    sample_col="Sample",
    category_col="cn",
)
```

Default figure formats:

| Function | Default |
|---|---|
| `plot_spatial()` | PNG |
| `plot_stacked_bar()` | PDF |
| `plot_cell_boundaries()` | PNG |
| `plot_annotation_polygons()` | PNG |
| `plot_marker_positivity_heatmap()` | PDF |
| `plot_marker_intensity_heatmap()` | PDF |
| `plot_cn_heatmap()` | PDF |

The same `save_png` / `save_pdf` controls are available for stacked bars and
cell-boundary plots. Heatmaps also retain their existing SVG and TIFF controls:

```python
qxy.plot_stacked_bar(adata, save_png=False, save_pdf=True)
qxy.plot_cell_boundaries(adata, save_png=True, save_pdf=False)

# PNG and PDF together.
qxy.plot_marker_positivity_heatmap(
    adata,
    save_png=True,
    save_pdf=True,
    save_svg=False,
    save_tiff=False,
)
qxy.plot_cn_heatmap(
    adata,
    save_png=True,
    save_pdf=True,
    save_svg=False,
    save_tiff=False,
)
```

Defaults produce one publication-appropriate figure file: PDF for stacked bars
and heatmaps, and PNG for spatial, cell-boundary, and annotation-polygon plots.
Annotation-polygon QC plots remain PNG-only by design.

Use a short image label column (`ImageID`) instead of the full QuPath `Image` name:

```python
qxy.plot_spatial(adata, sample_col="ImageID")
```

Reads from `adata.obsm["spatial"]` and `adata.obs[category_col]`. Spatial plots
default to bounding-box centering (`center_method="bbox"`) so each sample is
framed by its full X/Y extent rather than by an asymmetric cell-density median.
Colour palette cached in `adata.uns["qxycell"]["palettes"]`.

Plot cell boundary polygons instead of centroid dots when
`adata.obs["cell_polygon_wkt"]` is available:

```python
qxy.plot_cell_boundaries(
    adata,
    sample_col="Sample",
    samples=["sample_A", "sample_B"],
    label_celltypes="Tumor",
    auto_figsize=True,
    flip_y=True,
    save_pdf=False,
)
```

Boundary plots are more memory-intensive than dot plots because each cell is
drawn as a polygon. Cell labels are optional and only added for the cell type(s)
specified by `label_celltypes`.

Plot the original QuPath annotation polygons in one figure per image:

```python
# Reloads GeoJSON from adata.uns["qxycell"]["project_dir"] by default.
qxy.plot_annotation_polygons(adata, show=False)

# Override the original location when the project has moved.
qxy.plot_annotation_polygons(adata, project_dir="/path/to/qupath_project")
```

QXYCell stores annotation membership, rather than annotation geometry, in the
AnnData object. This function therefore reloads polygon geometry from the
original project folder and uses the stored run pixel size for coordinate scaling.
Cell locations are displayed as a low-resolution density underlay by default,
which remains efficient for images containing hundreds of thousands of cells.
Polygons are drawn as boundaries only by default (`fill=False`); pass
`fill=True` to add translucent polygon fills.
Use `cell_underlay=False` for polygon-only plots or adjust the density grid with
`underlay_bins=256`. These QC plots are saved as PNG files only.

Other annotation-polygon controls include `images`, `colors`,
`underlay_cmap`, `underlay_alpha`, `fill_alpha`, `boundary_linewidth`,
`flip_y`, `figsize`, and `dpi`.

## Stacked bar plots

Cell type or CN frequency per sample:

```python
# Cell type frequency per sample
qxy.plot_stacked_bar(adata)

# Averaged by group
qxy.plot_stacked_bar(adata, group_col="group")

# Selected cell types only
qxy.plot_stacked_bar(adata, celltypes=["T_cell", "Macrophage", "B_cell"])

# CN frequency
qxy.plot_stacked_bar(adata, category_col="cn")

# Journal widths: "single" = 90 mm (default), "double" = 180 mm
qxy.plot_stacked_bar(adata, width="double")

# Wider bars (default 15 mm per bar)
qxy.plot_stacked_bar(adata, bar_width_mm=20)
```

Reads from `adata.obs[category_col]` and `adata.obs[sample_col]`. Colour palette cached in `adata.uns["qxycell"]["palettes"]`.

## Heatmaps

Marker positivity and intensity are separate functions and save PDF by default:

```python
# Fraction positive per cell type × marker (batlow colormap, 0–1)
qxy.plot_marker_positivity_heatmap(adata)

# Z-score mean intensity (coolwarm, centred 0, ±3)
qxy.plot_marker_intensity_heatmap(adata)  # defaults to markers actually thresholded

# Journal single-column (90 mm) or double-column (180 mm) width
qxy.plot_marker_positivity_heatmap(adata, width="single")
qxy.plot_marker_intensity_heatmap(adata, width="double")

# Selected markers, no column clustering
qxy.plot_marker_intensity_heatmap(adata, markers=["CD45", "CD4", "CD8"], cluster_cols=False)

# Add a cell type colour strip on the left
qxy.plot_marker_positivity_heatmap(adata, row_strip=True)

# Group rows by CN instead of cell type
qxy.plot_marker_positivity_heatmap(adata, category_col="cn")
```

Positivity mode reads `adata.obs["<marker>_pos"]` columns. Intensity mode reads `adata.X` (z-scored per marker column).
The intensity matrix starts with the mean intensity for each category × marker,
then Z-scores each marker across categories. It is not a median or count
heatmap. Positivity is the positive-cell count divided by all cells in that
category. Both functions save one PDF plus the plotted matrix as CSV by
default. The legacy `plot_marker_heatmap(values=...)` entry point remains only
for older code.

When `markers=None`, the default marker heatmaps use the exact classifier source
measurement recorded during thresholding, including whether QuPath used `Mean`
or `Median`, and display its canonical compartment-aware marker name. Explicit
`markers=[...]` selections retain their existing behavior.

CN abundance heatmaps:

```python
# CN fraction per Sample, falling back to Image — columns sum to 1 (colorbar: f ↓)
qxy.plot_cn_heatmap(adata)

# Force one heatmap column per QuPath image
qxy.plot_cn_heatmap(adata, sample_col="Image")

# CN composition per CN — rows sum to 1 (colorbar: f →)
qxy.plot_cn_heatmap(adata, normalize="cn")

# Both normalisations in one call
qxy.plot_cn_heatmap(adata, normalize="both")

# Colour strip above samples grouped by condition
qxy.plot_cn_heatmap(adata, condition_col="group")

# Use a custom CN/category column
qxy.plot_cn_heatmap(adata, category_col="cn")
```

Reads from `adata.obs["cn"]` by default, or `adata.obs[category_col]` when
supplied. When `sample_col` is omitted, usable `Sample` labels are preferred
and QXYCell falls back to `Image`; an explicit `sample_col` is always honored.
Rows and columns are reordered by hierarchical clustering (no dendrogram).
Pass `cluster_rows=False` or `cluster_cols=False` to preserve input order.

## Colormaps

QXYCell resolves Crameri scientific colormap short names automatically. If `cmcrameri` is not installed a comparable matplotlib colormap is substituted silently.

| Alias | Type | Good for |
|---|---|---|
| `"batlow"` | Sequential | Positivity, CN abundance (default) |
| `"roma"` | Diverging | Z-score intensity (blue→white→red) |
| `"vik"` | Diverging | Z-score intensity (cooler tones) |

```python
qxy.plot_marker_intensity_heatmap(adata, cmap="roma")
qxy.plot_cn_heatmap(adata, cmap="vik")
```

Heatmap tiles are vector paths in PDF and SVG output, so they remain sharp
when enlarged in publication layouts. PNG and TIFF remain raster formats.

Any standard matplotlib colormap name also works (`cmap="viridis"` etc.).

## Colour consistency

`plot_spatial`, `plot_stacked_bar`, and the marker heatmap functions (with `row_strip=True`) share a per-category colour palette cached in `adata.uns["qxycell"]["palettes"]`. The same cell type or CN always gets the same colour across all plot types. Cell type plots use the glasbey palette; CN plots use `tab20` by default.

To regenerate a palette after adding new categories:

```python
adata.uns["qxycell"]["palettes"].pop("celltype")  # or "cn"
```

## TMA and CoreID

TMA core identity enters QXYCell through the QuPath cell measurement export.
Use this workflow:

1. Create and label the TMA grid in QuPath using the commands under **TMA**.
2. Confirm that the detected cells are associated with the intended cores.
3. Choose **Measure > Export measurements**, select the project images, set
   **Export type** to cells, and choose a tab separator.
4. Export `measurements.tsv` and confirm that the selected columns include the
   exact column name `TMA Core`.
5. Run `qxy.import_cells()`. QXYCell preserves `TMA Core` and
   automatically creates categorical `CoreID`.

QuPath documents the [TMA grid commands](https://qupath.readthedocs.io/en/stable/docs/reference/commands.html#tma)
and the [project measurement exporter](https://qupath.readthedocs.io/en/stable/docs/tutorials/exporting_measurements.html).

```python
adata = qxy.import_cells(project_dir)
adata.obs[["TMA Core", "CoreID"]].head()
adata.obs["CoreID"].value_counts()
```

The output is a categorical `adata.obs["CoreID"]`, and the summary is stored in
`adata.uns["qxycell_core_ids_from_measurements"]` and
`adata.uns["qxycell"]["measurement_core_assignment"]`.

If the measurement table has no `TMA Core` column, QXYCell reports
zero CoreIDs and does not add `adata.obs["CoreID"]`.

## Save and load

`qxy.save(adata)` writes the current AnnData object to a compressed `.h5ad`
file. This preserves the full analysis state: `adata.X`, `adata.obs`,
`adata.var`, `adata.obsm` such as `adata.obsm["spatial"]`, and `adata.uns`
metadata including QXYCell output paths, palettes, annotation summaries, CN
results, and other analysis state already stored on the object.

By default, QXYCell saves to the active run folder:

```text
qupath_project_run_YYMMDD_HHMM/h5ad/qxycell.h5ad
```

If `adata.uns["qxycell"]["h5ad_path"]` already exists, `qxy.save(adata)` updates
that same file. Pass `output_dir=` or `path=` to save somewhere else.

```python
# Save to the current run's output folder
qxy.save(adata)

# Save to a chosen output folder
qxy.save(adata, output_dir="my_qxy_output")

# Save to an exact file path
qxy.save(adata, path="my_analysis.h5ad")

# Reload the most recent save
adata = qxy.load_latest()

# Reload a specific file
adata = qxy.load("path/to/qxycell.h5ad")
```

## AnnData structure summary

| Location | Added by | Contents |
|---|---|---|
| `adata.obs["Image"]` | `qxy.import_cells()` | QuPath image name per cell |
| `adata.obs["Xµm"]`, `adata.obs["Yµm"]` | `qxy.import_cells()` | Cell centroid x/y coordinates in microns |
| `adata.obs["<marker>_pos"]` | explicit threshold stages | Boolean marker positivity columns |
| `adata.obs["annotation__<label>"]` | `qxy.add_annotations()` | Boolean annotation membership columns |
| `adata.obs["cell_polygon_wkt"]` | `qxy.add_annotations()` / `qxy.load_cell_polygons()` | Cell segmentation polygon geometry as WKT strings |
| `adata.obs["Sample"]` | `qxy.add_annotations()` / `qxy.assign_samples()` | Sample label from annotations with `Sample` in the label |
| `adata.obs["TMA Core"]` | `qxy.import_cells()` | QuPath measurement-table TMA core label |
| `adata.obs["CoreID"]` | `qxy.import_cells()` / `qxy.assign_core_ids_from_measurements()` | Categorical CoreID derived only from `TMA Core` |
| `adata.obs["celltype"]` | `qxy.celltype()` | Assigned cell type string |
| `adata.obs["cn"]` | `qxy.cn_kmeans()` | CN cluster label (int, then renamed to string) |
| `adata.obs[*metadata cols*]` | `qxy.add_metadata()` | Sample metadata broadcast to all cells |
| `adata.X` | `qxy.import_cells()` | Marker intensity matrix (cells × markers) |
| `adata.var` | `qxy.import_cells()` | Marker names and metadata |
| `adata.obsm["spatial"]` | `qxy.import_cells()` | Cell centroid x/y coordinates in microns |
| `adata.obsm["cn_profile"]` | `qxy.cn_knn()` | Per-cell local cell type composition (sums to 1) |
| `adata.uns["qxycell"]` | staged functions | Run metadata, output paths, stage status, and provenance |
| `adata.uns["qxycell_annotation_labels"]` | `qxy.add_annotations()` | Annotation class → column name map |
| `adata.uns["qxycell_thresholding"]` | explicit threshold stages | Threshold source and positivity-column summary |
| `adata.uns["qxycell_sample_annotations"]` | `qxy.assign_samples()` | Sample assignment summary |
| `adata.uns["qxycell_core_ids_from_measurements"]` | `qxy.assign_core_ids_from_measurements()` | Measurement-derived CoreID assignment summary |
| `adata.uns["qxycell_qc"]` | `qxy.qc()` | QC metrics per sample |
| `adata.uns["qxycell_sample_metadata"]` | `qxy.add_metadata()` | Metadata match summary |
| `adata.uns["qxycell_celltyping"]` | `qxy.celltype()` | Cell typing rule summary |
| `adata.uns["cn"]` | `qxy.cn_knn()` / `qxy.cn_kmeans()` | CN run parameters, cell type list, label map |

## Optional workflow shortcut

The staged API above is recommended for reviewable, rerunnable analyses.
`qxy.workflow()` remains available when a single convenience call is preferred:

```python
adata = qxy.workflow(
    "/path/to/qupath_project",
    threshold_file="/path/to/reviewed_thresholds.tsv",
    apply_thresholds=True,
    sample_metadata="sample_metadata.tsv",
    sample_col="Image",
    celltype_logic="celltype_logic.yaml",
)
```
