# QXYCell

<p align="center">
  <img src="assets/qxycell-icon.png" alt="QXYCell icon" width="180">
</p>

QXYCell converts QuPath single-cell exports into analysis-ready AnnData `.h5ad` objects for cell typing, neighbourhood analysis, spatial analysis, and plotting.

## Documentation

Full function reference: https://jsonmad.github.io/QXYCell/

## Installation

### conda (recommended)

Clone the repo, then create and activate the environment from the repo root:

```bash
git clone https://github.com/jsonmad/QXYCell.git
cd QXYCell
conda env create -f environment.yml
conda activate qxycell
```

To update an existing environment after pulling new changes:

```bash
cd QXYCell
git pull
conda env update -f environment.yml --prune
```

### pip

```bash
pip install qxycell
```

This installs all required dependencies: `anndata`, `colorcet`, `cmcrameri`, `geopandas`, `matplotlib`, `numpy`, `pandas`, `pyyaml`, `scikit-learn`, `scipy`, `seaborn`, and `shapely`.

## Quick start

```python
import qxycell as qxy

# Validate your QuPath export before running
report = qxy.check("/path/to/qupath_export")

# Import all cells into an AnnData object
adata = qxy.run("/path/to/qupath_export")

# Apply marker thresholds to create <marker>_pos columns
threshold_summary = qxy.threshold(adata, "/path/to/qupath_export")

# Apply cell type rules after thresholding
celltype_summary = qxy.celltype(adata, "celltype_logic.yaml")
```

Both functions write outputs to a timestamped folder: `qxy_outputs_YYMMDD-HHMM/`.

**`qxy.run()` populates the AnnData as follows:**

| Location | Contents |
|---|---|
| `adata.obs` | Per-cell metadata: `Image`, `Object ID`, `Xµm`, `Yµm`, optional `TMA Core` / `Parent`, automatic `CoreID` when measurement core metadata exists, `annotation__<label>` boolean columns, `Sample` when sample annotations exist, and `cell_polygon_wkt` when cell GeoJSON is available |
| `adata.X` | Marker intensity matrix (cells × markers) |
| `adata.var` | Marker names and metadata |
| `adata.obsm["spatial"]` | Cell centroid coordinates in microns, shape `(n_cells, 2)` |
| `adata.uns["qxycell"]` | Run metadata: output directory, timestamps, file paths |
| `adata.uns["qxycell_annotation_labels"]` | Map of annotation class names to `annotation__<label>` column names |

## QuPath inputs

QXYCell is built around manual QuPath exports. Required files:

- **Cell measurement table** — `measurements.csv` or `measurements.tsv` exported from QuPath. One table may contain cells from multiple images.
- **Threshold TSV/CSV** — the source of truth for marker positivity thresholds. Use a filled table such as `thresholds.tsv` or `thresholds_YYMMDD-HHMM.tsv`.
- **Object classifier JSONs** *(template source)* — single-measurement classifiers saved under `classifiers/object_classifiers/*.json`. QXYCell can convert these JSONs into a fresh timestamped threshold table, but existing threshold tables remain the active source.
- **Annotation GeoJSON** — exported QuPath annotation polygons, with measurements excluded. Regular annotation classification/name labels become boolean `annotation__<label>` columns in `adata.obs`. Annotations with `Sample` in the label define sample boundaries and are collapsed into one `adata.obs["Sample"]` column; annotations labelled `Ignore` mark regions to exclude. Annotation labels that exactly match measurement-derived core IDs are reported by `qxy.check()` but are not kept as annotation columns by `qxy.run()`, because measurement `CoreID` is preferred.
- **Cell segmentation GeoJSON** *(optional)* — exported cell objects for all cells, measurements excluded. Provides geometry for spatial analysis.
- **TMA core GeoJSON** *(optional)* — TMA core boundaries for TMA projects. These are not assigned by default; call `qxy.assign_tma_cores()` explicitly when geometry-based core assignment is needed.

Required measurement columns: `Image`, `Object ID`, `Centroid X µm`, `Centroid Y µm`.
These QuPath centroid columns are stored in `adata.obs` as `Xµm` and `Yµm`.
Optional measurement columns `TMA Core` and `Parent` are preserved when present
and are automatically collapsed into `adata.obs["CoreID"]` by `qxy.run()`.

Threshold tables are the source used by `qxy.check()` and `qxy.threshold()`.
`qxy.check()` always writes a fill-in threshold table based on the measurement
columns it finds:

```text
qxy_outputs_YYMMDD-HHMM/tables/thresholds_YYMMDD-HHMM.tsv
```

The table includes only measurement columns whose names contain `mean` or
`median`. Each row is one feature, and each image gets its own threshold
column.

When usable object classifier JSONs are present, `qxy.check()` writes a fresh
timestamped threshold table into the output `tables/` folder. This generated
table is a template/audit artifact; it does not overwrite any threshold table
in the QuPath export folder. To generate a fresh threshold table from JSONs
explicitly, call:

```python
threshold_path = qxy.generate_threshold_table("/path/to/qupath_export")
```

Object classifier JSONs never modify an existing threshold table. To manually
edit thresholds, fill in or edit the image threshold columns, save the finished
file back into the QuPath export folder, then rerun `qxy.check()` or
`qxy.threshold()`. By default, `qxy.threshold()` uses the newest recognized
threshold table in the QuPath export folder. To force a specific table, pass
`threshold_file="path/to/thresholds.tsv"` to `qxy.check()` or `qxy.threshold()`.
The compact `CheckReport` summary includes `Threshold source: ...`, and
`check_report.txt` includes both `Active threshold source: ...` and
`Generated threshold template: ...` so the source used for marker calls is
separate from the newly generated template.

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

## Annotations

Remove cells inside `Ignore` regions (tissue folds, artefacts):

```python
adata = qxy.remove_ignore(adata)
```

Removes rows from `adata.obs` in-place. Cells inside any `annotation__Ignore` polygon are dropped.

`qxy.run()` converts annotations with `Sample` in the label into one `Sample`
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

Results stored in `adata.uns["qxycell_qc"]`. HTML report and TSV tables written to `qxy_outputs_YYMMDD-HHMM/qc/`.

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

Generate a first-pass cell type prompt for an LLM:

```python
prompt = qxy.celltype_prompt(adata)
print(prompt)
```

This prints the prompt, saves it to `qxy_outputs_YYMMDD-HHMM/celltype/`, and returns the string. After editing and saving a YAML logic file in that folder, apply it:

```python
summary = qxy.celltype(adata)
```

Adds `adata.obs["celltype"]` (string). Cell type assignment summary stored in `adata.uns["qxycell_celltyping"]`. QXYCell applies the most recently saved YAML in the celltype folder and prints the file path it used.

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

Labels are derived from `adata.obsm["cn_profile"]` using a priority rule: clusters where one cell type dominates (≥ 50 %) are labelled `<type> hi`; mixed clusters show the top two contributors (e.g. `CD8 T/Macrophage`). Renames values in `adata.obs["cn"]` (string). Label map stored in `adata.uns["cn"]["label_map"]`. Summary saved to `qxy_outputs_YYMMDD-HHMM/cn/cn_labels.csv`.

Optionally shorten long cell type names before labelling:

```python
label_table = qxy.cn_name(adata, compaction={"CD8+PD1+LAG3+": "PD1 LAG3 CD8"})
```

## Spatial plots

Plot cell type distributions in tissue space:

```python
# One plot per sample
qxy.plot_spatial(adata)

# Combined multi-panel figure
qxy.plot_spatial(adata, combined=True, max_cols=3)

# Selected samples only
qxy.plot_spatial(adata, samples=["sample_A", "sample_B"])

# Fixed square window (microns)
qxy.plot_spatial(adata, fixed_window_um=11500)

# Figure panel aspect from selected sample X/Y extent
qxy.plot_spatial(adata, auto_figsize=True)

# Use raw coordinate orientation instead of the default image-viewer y flip
qxy.plot_spatial(adata, flip_y=False)

# Plot CNs instead of cell types
qxy.plot_spatial(adata, category_col="cn")
```

Use a short image label column (`ImageID`) instead of the full QuPath `Image` name:

```python
qxy.plot_spatial(adata, sample_col="ImageID")
```

Reads from `adata.obsm["spatial"]` and `adata.obs[category_col]`. Colour palette cached in `adata.uns["qxycell"]["palettes"]`.

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

# Wider bars (default 8 mm per bar)
qxy.plot_stacked_bar(adata, bar_width_mm=12)
```

Reads from `adata.obs[category_col]` and `adata.obs[sample_col]`. Colour palette cached in `adata.uns["qxycell"]["palettes"]`.

## Heatmaps

Marker positivity and intensity heatmaps — saves PDF, SVG, and TIFF (600 dpi):

```python
# Fraction positive per cell type × marker (batlow colormap, 0–1)
qxy.plot_marker_heatmap(adata)

# Z-score mean intensity (coolwarm, centred 0, ±3)
qxy.plot_marker_heatmap(adata, values="intensity")

# Both in one call
qxy.plot_marker_heatmap(adata, values="both")

# Journal single-column (90 mm) or double-column (180 mm) width
qxy.plot_marker_heatmap(adata, width="single")
qxy.plot_marker_heatmap(adata, width="double")

# Selected markers, no column clustering
qxy.plot_marker_heatmap(adata, markers=["CD45", "CD4", "CD8"], cluster_cols=False)

# Add a cell type colour strip on the left
qxy.plot_marker_heatmap(adata, row_strip=True)

# Group rows by CN instead of cell type
qxy.plot_marker_heatmap(adata, category_col="cn")
```

Positivity mode reads `adata.obs["<marker>_pos"]` columns. Intensity mode reads `adata.X` (z-scored per marker column).

CN abundance heatmaps:

```python
# CN fraction per sample — columns sum to 1 (colorbar: f ↓)
qxy.plot_cn_heatmap(adata)

# CN composition per CN — rows sum to 1 (colorbar: f →)
qxy.plot_cn_heatmap(adata, normalize="cn")

# Both normalisations in one call
qxy.plot_cn_heatmap(adata, normalize="both")

# Colour strip above samples grouped by condition
qxy.plot_cn_heatmap(adata, condition_col="group")

# Use a custom CN/category column
qxy.plot_cn_heatmap(adata, category_col="cn")
```

Reads from `adata.obs["cn"]` by default, or `adata.obs[category_col]` when supplied, plus `adata.obs[sample_col]`. Rows and columns are reordered by hierarchical clustering (no dendrogram). Pass `cluster_rows=False` or `cluster_cols=False` to preserve input order.

## Colormaps

QXYCell resolves Crameri scientific colormap short names automatically. If `cmcrameri` is not installed a comparable matplotlib colormap is substituted silently.

| Alias | Type | Good for |
|---|---|---|
| `"batlow"` | Sequential | Positivity, CN abundance (default) |
| `"roma"` | Diverging | Z-score intensity (blue→white→red) |
| `"vik"` | Diverging | Z-score intensity (cooler tones) |

```python
qxy.plot_marker_heatmap(adata, values="intensity", cmap="roma")
qxy.plot_cn_heatmap(adata, cmap="vik")
```

Any standard matplotlib colormap name also works (`cmap="viridis"` etc.).

## Colour consistency

`plot_spatial`, `plot_stacked_bar`, and `plot_marker_heatmap` (with `row_strip=True`) share a per-category colour palette cached in `adata.uns["qxycell"]["palettes"]`. The same cell type or CN always gets the same colour across all plot types. Cell type plots use the glasbey palette; CN plots use `tab20` by default.

To regenerate a palette after adding new categories:

```python
adata.uns["qxycell"]["palettes"].pop("celltype")  # or "cn"
```

## TMA and CoreID

For QuPath TMA dearrayer exports, use measurement-derived core IDs as the
default CoreID path. `qxy.run()` keeps the optional measurement columns
`"TMA Core"` and `"Parent"` when they are present and automatically collapses
them into one `adata.obs["CoreID"]` column.

```python
adata = qxy.run(project_dir)
adata.obs["CoreID"].value_counts()
```

By default this uses `"TMA Core"` first and then falls back to `"Parent"` for
cells where `"TMA Core"` is missing. The output is a categorical
`adata.obs["CoreID"]`, and the summary is stored in
`adata.uns["qxycell_core_ids_from_measurements"]` and
`adata.uns["qxycell"]["measurement_core_assignment"]`.

When annotation GeoJSON contains labels that exactly match these measurement
core IDs, `qxy.check()` reports them as "GeoJSON derived TMA CoreIDs".
`qxy.run()` then uses only the measurement `CoreID` values and skips those
matching GeoJSON labels so they do not appear as ordinary
`annotation__<core>` columns.

Use GeoJSON-based TMA assignment separately only when the core labels need to
come from core boundary geometry rather than from the measurement CSV:

```python
tma = qxy.assign_tma_cores(
    adata,
    "geojson/tma_cores.geojson",
    sample_col="Image",
)
```

Adds `adata.obs["tma_core"]` (string core label, `Unassigned` for unassigned
cells). Cells whose centroids fall inside overlapping TMA cores are left as
`Unassigned` rather than assigned arbitrarily. Summary stored in
`adata.uns["qxycell_tma"]`. GeoJSON files are matched to `adata.obs[sample_col]`
using the filename stem.

## Save and load

`qxy.save(adata)` writes the current AnnData object to a compressed `.h5ad`
file. This preserves the full analysis state: `adata.X`, `adata.obs`,
`adata.var`, `adata.obsm` such as `adata.obsm["spatial"]`, and `adata.uns`
metadata including QXYCell output paths, palettes, annotation summaries, CN
results, and other analysis state already stored on the object.

By default, QXYCell saves to the active run folder:

```text
qxy_outputs_YYMMDD-HHMM/run/h5ad/qxycell_YYMMDD-HHMM.h5ad
```

If `adata.uns["qxycell"]["h5ad_path"]` already exists, `qxy.save(adata)` updates
that same file. Pass `output_dir=` or `path=` to save somewhere else.

```python
# Save to the current run's output folder
qxy.save(adata)

# Save to a chosen output folder
qxy.save(adata, output_dir="qxy_outputs_manual")

# Save to an exact file path
qxy.save(adata, path="my_analysis.h5ad")

# Reload the most recent save
adata = qxy.load_latest()

# Reload a specific file
adata = qxy.load("path/to/qxycell_YYMMDD-HHMM.h5ad")
```

## AnnData structure summary

| Location | Added by | Contents |
|---|---|---|
| `adata.obs["Image"]` | `qxy.run()` | QuPath image name per cell |
| `adata.obs["Xµm"]`, `adata.obs["Yµm"]` | `qxy.run()` | Cell centroid x/y coordinates in microns |
| `adata.obs["<marker>_pos"]` | `qxy.threshold()` / `qxy.apply_thresholds()` | Boolean marker positivity columns |
| `adata.obs["annotation__<label>"]` | `qxy.run()` | Boolean annotation membership columns |
| `adata.obs["cell_polygon_wkt"]` | `qxy.run()` / `qxy.load_cell_polygons()` | Cell segmentation polygon geometry as WKT strings |
| `adata.obs["Sample"]` | `qxy.run()` / `qxy.assign_samples()` | Sample label from annotations with `Sample` in the label |
| `adata.obs["TMA Core"]`, `adata.obs["Parent"]` | `qxy.run()` | Optional QuPath measurement metadata columns when present |
| `adata.obs["CoreID"]` | `qxy.run()` / `qxy.assign_core_ids_from_measurements()` | Core ID from measurement metadata |
| `adata.obs["celltype"]` | `qxy.celltype()` | Assigned cell type string |
| `adata.obs["cn"]` | `qxy.cn_kmeans()` | CN cluster label (int, then renamed to string) |
| `adata.obs["tma_core"]` | `qxy.assign_tma_cores()` | TMA core label |
| `adata.obs[*metadata cols*]` | `qxy.add_metadata()` | Sample metadata broadcast to all cells |
| `adata.X` | `qxy.run()` | Marker intensity matrix (cells × markers) |
| `adata.var` | `qxy.run()` | Marker names and metadata |
| `adata.obsm["spatial"]` | `qxy.run()` | Cell centroid x/y coordinates in microns |
| `adata.obsm["cn_profile"]` | `qxy.cn_knn()` | Per-cell local cell type composition (sums to 1) |
| `adata.uns["qxycell"]` | `qxy.run()` | Run metadata, output paths, colour palettes |
| `adata.uns["qxycell_annotation_labels"]` | `qxy.run()` | Annotation class → column name map |
| `adata.uns["qxycell_thresholding"]` | `qxy.threshold()` / `qxy.apply_thresholds()` | Threshold source and positivity-column summary |
| `adata.uns["qxycell_sample_annotations"]` | `qxy.assign_samples()` | Sample assignment summary |
| `adata.uns["qxycell_core_ids_from_measurements"]` | `qxy.assign_core_ids_from_measurements()` | Measurement-derived CoreID assignment summary |
| `adata.uns["qxycell_qc"]` | `qxy.qc()` | QC metrics per sample |
| `adata.uns["qxycell_sample_metadata"]` | `qxy.add_metadata()` | Metadata match summary |
| `adata.uns["qxycell_celltyping"]` | `qxy.celltype()` | Cell typing rule summary |
| `adata.uns["cn"]` | `qxy.cn_knn()` / `qxy.cn_kmeans()` | CN run parameters, cell type list, label map |
| `adata.uns["qxycell_tma"]` | `qxy.assign_tma_cores()` | TMA assignment summary |

## Workflow shortcut

Run the full pipeline in one call:

```python
adata = qxy.workflow(
    "/path/to/qupath_export",
    sample_metadata="sample_metadata.tsv",
    sample_col="Image",
    celltype_logic="celltype_logic.yaml",
)
```
