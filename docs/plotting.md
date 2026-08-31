# Plotting

QXYCell provides spatial, categorical, and heatmap functions for cell types,
cellular neighbourhoods, marker measurements, annotations, and cell
boundaries. Plotting functions reuse the active output folder and consistent
category palettes stored in AnnData.

## Spatial plots

```python
# Prefer Sample when available; otherwise use Image.
qxy.plot_spatial(adata)

# Force one plot per QuPath image or Sample annotation.
qxy.plot_spatial(adata, sample_col="Image")
qxy.plot_spatial(adata, sample_col="Sample")

# Plot cellular neighbourhoods instead of cell types.
qxy.plot_spatial(adata, category_col="cn")
```

When `sample_col` is omitted, usable `Sample` labels are preferred and QXYCell
falls back to `Image`. Pass a value explicitly to force the grouping.
`samples=` always selects values from the resolved sample column.

Cells with missing sample values are excluded by default. Use
`include_missing_samples=True` only when a separate `nan` panel is wanted.

Common layout controls:

```python
qxy.plot_spatial(adata, combined=True, max_cols=3)
qxy.plot_spatial(adata, samples=["image_A.ome.tif", "image_B.ome.tif"])
qxy.plot_spatial(adata, fixed_window_um=11500)
qxy.plot_spatial(adata, center_method="median")
qxy.plot_spatial(adata, auto_figsize=True)
qxy.plot_spatial(adata, flip_y=False)
```

Spatial plots default to bounding-box centering and image-viewer y-axis
orientation. Coordinates come from `adata.obsm["spatial"]`.

Use a complete AnnData object as a grey underlay while plotting categories
from a filtered object:

```python
qxy.plot_spatial(
    adata_cn,
    underlay_adata=adata,
    sample_col="Sample",
    category_col="cn",
)
```

## Cell boundaries

When `adata.obs["cell_polygon_wkt"]` is available, plot cell polygons instead
of centroid dots:

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

Boundary plots are more memory-intensive because every cell is drawn as a
polygon. Labels are optional and are added only for the cell types selected by
`label_celltypes`.

## Annotation polygons

Plot the original QuPath annotation polygons:

```python
# Reload GeoJSON from the stored project folder.
qxy.plot_annotation_polygons(adata, show=False)

# Override the source when the project has moved.
qxy.plot_annotation_polygons(
    adata,
    project_dir="/path/to/qupath_project",
)
```

QXYCell stores annotation membership rather than polygon geometry in AnnData,
so this function reloads GeoJSON from the project folder and applies the saved
pixel size. A low-resolution cell-density underlay is enabled by default.
Polygons use boundary-only rendering unless `fill=True`.

Useful controls include `images`, `colors`, `underlay_cmap`,
`underlay_alpha`, `fill_alpha`, `boundary_linewidth`, `cell_underlay`,
`underlay_bins`, `flip_y`, `figsize`, and `dpi`. Annotation-polygon diagnostic
plots are PNG-only.

## Figure formats

Default outputs are:

| Function | Default format |
|---|---|
| `plot_spatial()` | PNG |
| `plot_stacked_bar()` | PDF |
| `plot_cell_boundaries()` | PNG |
| `plot_annotation_polygons()` | PNG |
| `plot_marker_positivity_heatmap()` | PDF |
| `plot_marker_intensity_heatmap()` | PDF |
| `plot_cn_heatmap()` | PDF |

Most plotting functions accept `save_png` and `save_pdf`:

```python
qxy.plot_spatial(adata, save_png=True, save_pdf=False)
qxy.plot_stacked_bar(adata, save_png=False, save_pdf=True)
qxy.plot_cell_boundaries(adata, save_png=True, save_pdf=True)
```

Heatmaps also support SVG and TIFF:

```python
qxy.plot_marker_positivity_heatmap(
    adata,
    save_png=True,
    save_pdf=True,
    save_svg=False,
    save_tiff=False,
)
```

## Stacked bar plots

Plot cell-type or neighbourhood frequencies per sample:

```python
qxy.plot_stacked_bar(adata)
qxy.plot_stacked_bar(adata, group_col="group")
qxy.plot_stacked_bar(
    adata,
    celltypes=["T_cell", "Macrophage", "B_cell"],
)
qxy.plot_stacked_bar(adata, category_col="cn")
```

Publication-width controls use `width="single"` for 90 mm or
`width="double"` for 180 mm. Set `bar_width_mm` to control individual bar
width.

## Marker heatmaps

Marker positivity and intensity are separate functions:

```python
# Fraction positive per category and marker, scaled 0–1.
qxy.plot_marker_positivity_heatmap(adata)

# Z-scored mean intensity per category and marker.
qxy.plot_marker_intensity_heatmap(adata)
```

Positivity reads `adata.obs["<marker>_pos"]`. Intensity starts with the mean
intensity for each category and marker, then Z-scores each marker across
categories. It is not a median or count heatmap.

```python
qxy.plot_marker_intensity_heatmap(
    adata,
    markers=["CD45", "CD4", "CD8"],
    cluster_cols=False,
)
qxy.plot_marker_positivity_heatmap(adata, row_strip=True)
qxy.plot_marker_positivity_heatmap(adata, category_col="cn")
```

When `markers=None`, heatmaps use the exact threshold source measurement,
including Mean or Median, and display the canonical compartment-aware marker
name. Each heatmap saves the plotted matrix as CSV alongside the figure. The
legacy `plot_marker_heatmap(values=...)` entry point remains for older code.

## Cellular-neighbourhood heatmaps

```python
# CN fraction per Sample, falling back to Image; columns sum to 1.
qxy.plot_cn_heatmap(adata)

# Force one column per QuPath image.
qxy.plot_cn_heatmap(adata, sample_col="Image")

# Cell-type composition per CN; rows sum to 1.
qxy.plot_cn_heatmap(adata, normalize="cn")

# Produce both normalizations.
qxy.plot_cn_heatmap(adata, normalize="both")

# Add a condition colour strip.
qxy.plot_cn_heatmap(adata, condition_col="group")
```

Rows and columns are reordered by hierarchical clustering without a
dendrogram. Use `cluster_rows=False` or `cluster_cols=False` to preserve input
order.

## Colormaps

QXYCell resolves Crameri scientific colormap aliases automatically. If
`cmcrameri` is unavailable, a comparable Matplotlib colormap is substituted.

| Alias | Type | Typical use |
|---|---|---|
| `batlow` | Sequential | Positivity and CN abundance |
| `roma` | Diverging | Z-scored intensity |
| `vik` | Diverging | Z-scored intensity with cooler tones |

```python
qxy.plot_marker_intensity_heatmap(adata, cmap="roma")
qxy.plot_cn_heatmap(adata, cmap="vik")
```

Standard Matplotlib names such as `viridis` also work. Heatmap tiles remain
vector paths in PDF and SVG output; PNG and TIFF are raster formats.

## Colour consistency

Spatial plots, stacked bars, and marker heatmaps with `row_strip=True` share a
per-category palette stored in `adata.uns["qxycell"]["palettes"]`. The same
cell type or neighbourhood keeps the same colour across plot types. Cell types
use a glasbey palette; neighbourhoods use `tab20` by default.

Regenerate a palette after adding categories:

```python
adata.uns["qxycell"]["palettes"].pop("celltype")  # or "cn"
```

## Related documentation

- Add plot grouping fields with the [sample-metadata guide](metadata.md).
- Prepare reviewed labels with the [cell-typing guide](cell_typing.md).
- Build neighbourhood labels with the
  [cellular-neighbourhood guide](cellular_neighbourhoods.md).
- Review the stored object in [AnnData and outputs](anndata_and_outputs.md).
- See every plotting parameter in the
  [function reference](QXYCell_function_reference.html).
