# QUXYCell

QUXYCell is a Python package for processing manually exported QuPath single-cell projects.

The first user-facing workflow is Python-first:

```python
import quxycell as qxy

report = qxy.check("/path/to/qupath_export")
adata = qxy.run("/path/to/qupath_export")
```

Version 1 assumes the user has already exported files from QuPath:

- `measurements.csv` or `measurements.tsv`
- `classifiers/object_classifiers/*.json`
- exported annotation GeoJSON files

Required measurement columns:

- `Image`
- `Object ID`
- `Centroid X µm`
- `Centroid Y µm`

QUXYCell does not parse QuPath `.qpdata` directly in v1.

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

QUXYCell prints the YAML file path it used.
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
