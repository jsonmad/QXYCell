# QuXYCell

<p align="center">
  <img src="assets/quxycell-icon.png" alt="QuXYCell icon" width="180">
</p>

QuXYCell is a Python package for processing manually exported QuPath single-cell projects.

## Documentation

The GitHub Pages documentation site is available at:

https://jsonmad.github.io/QuXYCell/

The first user-facing workflow is Python-first:

```python
import quxycell as qxy

report = qxy.check("/path/to/qupath_export")
adata = qxy.run("/path/to/qupath_export")
```

## QuPath Inputs

QuXYCell v1 is built around manual QuPath exports. It does not parse QuPath
`.qpdata` files directly.

Required QuPath products:

- Cell measurement table: one or more `measurements.csv` or `measurements.tsv`
  files exported from QuPath. A single table may contain cells from multiple
  images.
- Object classifier JSON files: QuPath classifier definitions, usually under
  `classifiers/object_classifiers/*.json`. QuXYCell reads threshold rules from
  these files and creates marker positivity columns in `adata.obs` named
  `<marker>_pos`.
- Annotation GeoJSON files: exported QuPath annotations containing polygon
  boundaries and annotation names/classes. QuXYCell assigns cells inside these
  polygons to annotation columns in `adata.obs`.

Required measurement columns:

- `Image`: the QuPath image name. QuXYCell uses this to keep cells associated
  with the correct image/sample.
- `Object ID`: the QuPath cell/object identifier.
- `Centroid X µm`: cell centroid x-coordinate in microns.
- `Centroid Y µm`: cell centroid y-coordinate in microns.

Optional QuPath products:

- Cell segmentation GeoJSON files: if exported, these can be used to preserve
  cell/object geometries where available.
- TMA core GeoJSON files: if exported, QuXYCell can assign cells to
  non-overlapping TMA core boundaries and preserve core metadata.

QuXYCell assumes the standard QuPath `Image` column name. Alternative aliases
are not used in v1.

If `output_dir` is omitted:

- `qxy.check(...)` writes reports to `qxy_outputs_YYMMDD-HHMM/`
- `qxy.run(...)` writes outputs to `qxy_outputs_YYMMDD-HHMM/`

To draft a first-pass cell type YAML prompt from a loaded AnnData object:

```python
prompt = qxy.celltype_prompt(adata)
print(prompt)
```

By default, this prints the prompt below the notebook cell, prints usage
instructions, and saves a timestamped prompt to the current run's
`qxy_outputs_YYMMDD-HHMM/celltype/` folder.
Use `save=False` to return the prompt without writing a file, or
`print_prompt=False` to avoid printing.

To apply cell typing from the newest YAML in the current run's
`qxy_outputs_YYMMDD-HHMM/celltype/` folder:

```python
summary = qxy.celltype(adata)
```

QuXYCell prints the YAML file path it used.
It also saves a cell type rule summary table to
`qxy_outputs_YYMMDD-HHMM/celltype/celltype_rules_summary_<logic-name>.tsv`.

Generic plotting helpers:

```python
qxy.plot_stacked_bar(adata, category_col="celltype", group_col="Group")
qxy.plot_spatial_celltypes(adata, category_col="celltype")
```

Use another image label column, such as a shortened `ImageID`, with
`sample_col`:

```python
qxy.plot_stacked_bar(adata, category_col="celltype", sample_col="ImageID")
qxy.plot_spatial_celltypes(adata, category_col="celltype", sample_col="ImageID")
```

Spatial plots use a shared centered x/y range across all selected samples,
centered by median cell coordinate by default, and include a 1 mm scale bar.
To force a fixed square window:

```python
qxy.plot_spatial_celltypes(adata, sample_col="ImageID", fixed_window_um=11500)
```

Remove cells inside annotation columns containing `Ignore`:

```python
adata = qxy.remove_ignore(adata)
```

Save and reload:

```python
qxy.save(adata)
adata = qxy.load_latest()
adata = qxy.load("qxy_outputs_YYMMDD-HHMM")
```

Add sample metadata:

```python
qxy.add_metadata(
    adata,
    "sample_metadata.tsv",
    sample_col="ImageID",
)
```

`sample_col` can be any column in `adata.obs`; the default is QuPath `Image`.
If the metadata table uses a different key column name, set
`metadata_sample_col`.

Generate QC tables and an HTML report:

```python
qc = qxy.qc(adata, sample_col="ImageID")
```

Assign non-overlapping TMA cores from GeoJSON:

```python
tma = qxy.assign_tma_cores(
    adata,
    "geojson/tma_cores.geojson",
    sample_col="Image",
)
```

TMA GeoJSON files are matched to `adata.obs[sample_col]` using the GeoJSON
filename stem. The default `sample_col="Image"` matches QuPath image names.

Use plot controls:

```python
qxy.plot_stacked_bar(
    adata,
    sample_col="ImageID",
    celltypes=["T_cell", "Macrophage"],
    palette="tab20",
)

qxy.plot_spatial_celltypes(
    adata,
    sample_col="ImageID",
    celltypes=["T_cell", "Macrophage"],
    samples=["beat-2", "beat-5"],
    combined=True,
    max_cols=2,
    palette="tab20",
)
```

Run the common notebook workflow in one command:

```python
adata = qxy.workflow(
    "/path/to/qupath_export",
    sample_metadata="sample_metadata.tsv",
    sample_col="Image",
    celltype_logic="celltype_logic.yaml",
)
```
