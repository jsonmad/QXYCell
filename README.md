# QuXYCell

<p align="center">
  <img src="assets/quxycell-icon.png" alt="QuXYCell icon" width="180">
</p>

QuXYCell converts manually exported QuPath single-cell projects into analysis-ready AnnData `.h5ad` objects, with built-in tools for cell typing, neighbourhood analysis, and publication-quality plotting.

## Documentation

Full function reference: https://jsonmad.github.io/QuXYCell/

## Installation

```bash
pip install quxycell
```

This installs all required dependencies: `anndata`, `cmcrameri`, `geopandas`, `matplotlib`, `numpy`, `pandas`, `pyyaml`, `scikit-learn`, `scipy`, `seaborn`, and `shapely`.

## Quick start

```python
import quxycell as qxy

# Validate your QuPath export before running
report = qxy.check("/path/to/qupath_export")

# Import all cells into an AnnData object
adata = qxy.run("/path/to/qupath_export")
```

Both functions write outputs to a timestamped folder: `qxy_outputs_YYMMDD-HHMM/`.

**`qxy.run()` populates the AnnData as follows:**

| Location | Contents |
|---|---|
| `adata.obs` | Per-cell metadata: `Image`, `Object ID`, `Centroid X µm`, `Centroid Y µm`, `<marker>_pos` boolean columns, `annotation__<label>` boolean columns |
| `adata.X` | Marker intensity matrix (cells × markers) |
| `adata.var` | Marker names and metadata |
| `adata.obsm["spatial"]` | Cell centroid coordinates in microns, shape `(n_cells, 2)` |
| `adata.uns["quxycell"]` | Run metadata: output directory, timestamps, file paths |
| `adata.uns["quxycell_annotation_labels"]` | Map of annotation class names to `annotation__<label>` column names |

## QuPath inputs

QuXYCell is built around manual QuPath exports. Required files:

- **Cell measurement table** — `measurements.csv` or `measurements.tsv` exported from QuPath. One table may contain cells from multiple images.
- **Object classifier JSONs** — single-measurement classifiers saved under `classifiers/object_classifiers/*.json`. QuXYCell reads the marker name and positivity threshold from each JSON and creates `<marker>_pos` boolean columns in `adata.obs`.
- **Annotation GeoJSON** — exported QuPath annotation polygons, with measurements excluded. Each annotation's class/name becomes a boolean `annotation__<label>` column in `adata.obs`. Annotations labelled `Sample-` define sample boundaries within an image; annotations labelled `Ignore` mark regions to exclude.
- **Cell segmentation GeoJSON** *(optional)* — exported cell objects for all cells, measurements excluded. Provides geometry for spatial analysis.
- **TMA core GeoJSON** *(optional)* — TMA core boundaries for TMA projects.

Required measurement columns: `Image`, `Object ID`, `Centroid X µm`, `Centroid Y µm`.

## Annotations

Remove cells inside `Ignore` regions (tissue folds, artefacts):

```python
adata = qxy.remove_ignore(adata)
```

Removes rows from `adata.obs` in-place. Cells inside any `annotation__Ignore` polygon are dropped.

Convert `Sample-` annotations into a `Sample` column:

```python
sample_summary = qxy.assign_samples(adata)
adata.obs["Sample"].value_counts()
```

Adds `adata.obs["Sample"]` (string). Summary stored in `adata.uns["quxycell_sample_annotations"]`.

## QC

Generate per-sample QC tables and an HTML report:

```python
qc = qxy.qc(adata, sample_col="Image")
```

Results stored in `adata.uns["quxycell_qc"]`. HTML report and TSV tables written to `qxy_outputs_YYMMDD-HHMM/qc/`.

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

Each column in the metadata file is broadcast to every cell belonging to that sample and added as a new column in `adata.obs`. For example, a `group` column in the TSV becomes `adata.obs["group"]`. Summary stored in `adata.uns["quxycell_sample_metadata"]`.

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

Adds `adata.obs["celltype"]` (string). Cell type assignment summary stored in `adata.uns["quxycell_celltyping"]`. QuXYCell applies the most recently saved YAML in the celltype folder and prints the file path it used.

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

Assign human-readable names to CN clusters:

```python
qxy.cn_name(adata, {0: "Immune-rich", 1: "Stromal", 2: "Tumour core"})
```

Renames values in `adata.obs["cn"]` (string). Label map stored in `adata.uns["cn"]["label_map"]`.

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

# Plot CNs instead of cell types
qxy.plot_spatial(adata, category_col="cn")
```

Use a short image label column (`ImageID`) instead of the full QuPath `Image` name:

```python
qxy.plot_spatial(adata, sample_col="ImageID")
```

Reads from `adata.obsm["spatial"]` and `adata.obs[category_col]`. Colour palette cached in `adata.uns["quxycell"]["palettes"]`.

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

Reads from `adata.obs[category_col]` and `adata.obs[sample_col]`. Colour palette cached in `adata.uns["quxycell"]["palettes"]`.

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
```

Reads from `adata.obs["cn"]` and `adata.obs[sample_col]`. Rows and columns are reordered by hierarchical clustering (no dendrogram). Pass `cluster_rows=False` or `cluster_cols=False` to preserve input order.

## Colormaps

QuXYCell resolves Crameri scientific colormap short names automatically. If `cmcrameri` is not installed a comparable matplotlib colormap is substituted silently.

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

`plot_spatial`, `plot_stacked_bar`, and `plot_marker_heatmap` (with `row_strip=True`) share a per-category colour palette cached in `adata.uns["quxycell"]["palettes"]`. The same cell type or CN always gets the same colour across all plot types. Cell type plots use the glasbey palette; CN plots use `tab20` by default.

To regenerate a palette after adding new categories:

```python
adata.uns["quxycell"]["palettes"].pop("celltype")  # or "cn"
```

## TMA

Assign cells to TMA cores from a GeoJSON boundary file:

```python
tma = qxy.assign_tma_cores(
    adata,
    "geojson/tma_cores.geojson",
    sample_col="Image",
)
```

Adds `adata.obs["tma_core"]` (string core label, `NaN` for unassigned cells). Summary stored in `adata.uns["quxycell_tma"]`. GeoJSON files are matched to `adata.obs[sample_col]` using the filename stem.

## Save and load

```python
# Save to the current run's output folder
qxy.save(adata)

# Reload the most recent save
adata = qxy.load_latest()

# Reload a specific file
adata = qxy.load("path/to/quxycell_YYMMDD-HHMM.h5ad")
```

## AnnData structure summary

| Location | Added by | Contents |
|---|---|---|
| `adata.obs["Image"]` | `qxy.run()` | QuPath image name per cell |
| `adata.obs["<marker>_pos"]` | `qxy.run()` | Boolean marker positivity columns |
| `adata.obs["annotation__<label>"]` | `qxy.run()` | Boolean annotation membership columns |
| `adata.obs["Sample"]` | `qxy.assign_samples()` | Sample label from `Sample-` annotations |
| `adata.obs["celltype"]` | `qxy.celltype()` | Assigned cell type string |
| `adata.obs["cn"]` | `qxy.cn_kmeans()` | CN cluster label (int, then renamed to string) |
| `adata.obs["tma_core"]` | `qxy.assign_tma_cores()` | TMA core label |
| `adata.obs[*metadata cols*]` | `qxy.add_metadata()` | Sample metadata broadcast to all cells |
| `adata.X` | `qxy.run()` | Marker intensity matrix (cells × markers) |
| `adata.var` | `qxy.run()` | Marker names and metadata |
| `adata.obsm["spatial"]` | `qxy.run()` | Cell centroid x/y coordinates in microns |
| `adata.obsm["cn_profile"]` | `qxy.cn_knn()` | Per-cell local cell type composition (sums to 1) |
| `adata.uns["quxycell"]` | `qxy.run()` | Run metadata, output paths, colour palettes |
| `adata.uns["quxycell_annotation_labels"]` | `qxy.run()` | Annotation class → column name map |
| `adata.uns["quxycell_sample_annotations"]` | `qxy.assign_samples()` | Sample assignment summary |
| `adata.uns["quxycell_qc"]` | `qxy.qc()` | QC metrics per sample |
| `adata.uns["quxycell_sample_metadata"]` | `qxy.add_metadata()` | Metadata match summary |
| `adata.uns["quxycell_celltyping"]` | `qxy.celltype()` | Cell typing rule summary |
| `adata.uns["cn"]` | `qxy.cn_knn()` / `qxy.cn_kmeans()` | CN run parameters, cell type list, label map |
| `adata.uns["quxycell_tma"]` | `qxy.assign_tma_cores()` | TMA assignment summary |

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
