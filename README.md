# QUXYCell

QUXYCell is a Python package for processing manually exported QuPath single-cell projects.

The first user-facing workflow is Python-first:

```python
import quxycell as qxy

report = qxy.check("/path/to/qupath_export", output_dir="./quxycell_output")
adata = qxy.run("/path/to/qupath_export", output_dir="./quxycell_output")
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

