# Analysis workflows

After importing measurements, annotations, and thresholds, QXYCell can add
descriptive summaries, sample metadata, reviewed cell-type assignments, and
cellular neighbourhoods to the active AnnData object.

## Dataset summary

Generate descriptive tables and an HTML report for the cells, markers,
annotations, samples, and cell types already present in AnnData:

```python
summary = qxy.dataset_summary(adata, sample_col="Sample")
```

Results are stored in `adata.uns["qxycell_dataset_summary"]` and written under
`dataset_summary/` in the active output folder, including:

- `dataset_summary.html`
- `dataset_overview.tsv`
- additional TSV tables for available dataset components

This is descriptive reporting. It does not validate image quality,
segmentation, staining, thresholds, batch effects, or spatial alignment.

## Sample metadata

Add sample-level metadata matched on image name:

```python
qxy.add_metadata(adata, "sample_metadata.tsv")
```

Use different matching columns when required:

```python
qxy.add_metadata(
    adata,
    "sample_metadata.tsv",
    sample_col="ImageID",          # column in adata.obs
    metadata_sample_col="sample",  # column in the TSV
)
```

Import selected columns only:

```python
qxy.add_metadata(
    adata,
    "sample_metadata.tsv",
    columns=["group", "mouse_id"],
)
```

Each metadata value is broadcast to every cell belonging to the matched
sample. The summary is stored in
`adata.uns["qxycell_sample_metadata"]`.

## Cell typing

Cell typing uses marker-positive calls rather than raw measurement values.
Apply thresholds first with exactly one explicit source:

```python
qxy.threshold_from_classifiers(adata)
# or
qxy.threshold_from_table(adata, "/path/to/thresholds.tsv")
```

For each marker, QXYCell compares the selected QuPath measurement with the
applicable global or per-image threshold. Values greater than or equal to the
threshold become `True` in `adata.obs["<MARKER>_pos"]`.

### Cell-type rule behavior

Cell-type YAML rules refer to positivity columns by marker name:

- `positive` requires every listed marker to be positive.
- `negative` requires every listed marker to be negative.
- `any_positive` requires at least one listed marker to be positive.
- Rules are evaluated from top to bottom; the first match wins.

Put specific cell types before broader populations so broad rules do not
capture them first. Compartment-aware classifier names such as `CD3-nuc` must
be used exactly in the YAML.

Generate the prompt used to draft a first-pass logic file:

```python
prompt = qxy.celltype_prompt(
    adata,
    context="Describe the tissue and expected populations",
)
print(prompt)
```

The prompt is returned, printed, and saved as
`celltype/current_prompt.txt`. Review the proposed YAML with a domain expert,
save it, and then apply it:

```python
summary = qxy.celltype(adata, "/path/to/celltype_logic.yaml")
```

This adds `adata.obs["celltype"]` and stores the rule summary in
`adata.uns["qxycell_celltyping"]`. If no path is supplied, QXYCell uses the
most recently saved YAML in the active `celltype/` folder and prints the path.

Rerunning thresholding removes threshold-dependent cell-type and feature
columns. Rerun `qxy.celltype()` afterwards. Rerunning cell typing replaces
columns owned by the previous YAML, including obsolete derived features.

## Cellular neighbourhoods

### Build local composition profiles

Use k-nearest neighbours to calculate the local cell-type composition around
each cell:

```python
qxy.cn_knn(adata, k=15)
```

This adds `adata.obsm["cn_profile"]`, a float32 matrix of shape
`(n_cells, n_cell_types)`. Parameters are stored in `adata.uns["cn"]`.

### Cluster neighbourhood profiles

```python
qxy.cn_kmeans(adata, n_cn=10)
```

This adds integer cluster labels in `adata.obs["cn"]` and updates
`adata.uns["cn"]` with clustering parameters.

### Assign descriptive names

```python
label_table = qxy.cn_name(adata)
```

Clusters dominated by one cell type (at least 50%) are labelled
`<type> hi`; mixed clusters show their top two contributors. The label map is
stored in `adata.uns["cn"]["label_map"]` and written to
`cn/cn_labels.csv`.

Long names can be shortened before labelling:

```python
label_table = qxy.cn_name(
    adata,
    compaction={"CD8+PD1+LAG3+": "PD1 LAG3 CD8"},
)
```

Generated labels are Windows-safe: path separators become ` + `, reserved
filename characters are replaced, repeated `..` segments are collapsed, and
trailing spaces or periods are removed.

## Next steps

- Create spatial and summary figures with the [plotting guide](plotting.md).
- Review stored fields and save/reload behavior in
  [AnnData and outputs](anndata_and_outputs.md).
- See all parameters in the
  [function reference](QXYCell_function_reference.html).
