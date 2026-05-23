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

- `qxy.check(...)` writes reports to `outputs/qxy_check/`
- `qxy.run(...)` writes outputs to `outputs/qxy_run/`

To draft a first-pass cell type YAML prompt from a loaded AnnData object:

```python
prompt = qxy.celltype_prompt(adata)
print(prompt)
```

By default, this prints the prompt below the notebook cell, prints usage
instructions, and saves a timestamped prompt to `outputs/qxy_run/celltype/`.
Use `save=False` to return the prompt without writing a file, or
`print_prompt=False` to avoid printing.

To apply cell typing from the newest YAML in `outputs/qxy_run/celltype/`:

```python
summary = qxy.celltype(adata)
```

QUXYCell prints the YAML file path it used.
