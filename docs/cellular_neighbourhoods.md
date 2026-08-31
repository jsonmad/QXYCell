# Cellular neighbourhoods

Cellular neighbourhood analysis describes the local mixture of cell types
around every cell. QXYCell builds a spatial k-nearest-neighbour composition
profile, clusters those profiles, and can assign deterministic descriptive
labels to the resulting neighbourhoods.

## Prerequisites

The AnnData object must contain:

- reviewed cell labels in `adata.obs["celltype"]`;
- cell coordinates in `adata.obsm["spatial"]`; and
- an image identifier in `adata.obs["Image"]`.

Run cell typing and inspect the assignments before neighbourhood analysis:

```python
adata.obs["celltype"].value_counts(dropna=False)
```

QXYCell builds neighbours independently within each image so cells from
different tissue sections cannot become neighbours. Keep `image_col="Image"`
unless another column also defines physically independent coordinate spaces.

## Decide which cells to include

By default, every cell type—including `Unknown`—contributes to local
composition profiles. Use the active object when all cells should be retained:

```python
adata_cn = adata
```

To exclude a category, subset a copy before both the KNN and clustering steps.
Give a filtered copy its own output folder immediately so it cannot replace the
full dataset checkpoint:

```python
adata_cn = adata[adata.obs["celltype"] != "Unknown"].copy()
qxy.save(
    adata_cn,
    output_dir="/path/to/qxycell_outputs/cn_without_unknown",
)
```

This changes the analysed cell population and therefore the neighbour graph.
Record and report any filtering decision.

## Build local composition profiles

For each cell, calculate the fraction of its nearest neighbours belonging to
each cell type:

```python
qxy.cn_knn(
    adata_cn,
    k=12,
    cell_type_col="celltype",
    image_col="Image",
)
```

The focal cell is excluded. If an image contains fewer than `k` other cells,
QXYCell uses the available neighbours. Results are stored as a float32 matrix
in `adata_cn.obsm["cn_profile"]`; each row contains the local cell-type
fractions for one cell. Parameters and the ordered cell-type list are stored in
`adata_cn.uns["cn"]`.

`k` controls the spatial scale of the analysis. Smaller values emphasize very
local organization and are more sensitive to individual cells; larger values
describe broader tissue context. Choose it according to cell density,
biological scale, and a documented sensitivity analysis.

## Cluster neighbourhood profiles

Cluster the local profiles with MiniBatchKMeans:

```python
qxy.cn_kmeans(
    adata_cn,
    n_cn=12,
    random_state=0,
)
```

This adds categorical cluster identifiers such as `N0`, `N1`, and `N2` to
`adata_cn.obs["cn"]`. The default random seed makes repeated runs with the same
data and parameters reproducible. `n_cn` controls the number of neighbourhood
classes and should be reviewed for stability and biological interpretability.

Use separate AnnData copies when comparing parameter choices so one result
does not overwrite another:

```python
adata_cn_8 = adata_cn.copy()
qxy.cn_kmeans(adata_cn_8, n_cn=8, random_state=0)

adata_cn_12 = adata_cn.copy()
qxy.cn_kmeans(adata_cn_12, n_cn=12, random_state=0)
```

## Assign descriptive names

Generate deterministic labels from the mean cell-type composition of each
cluster:

```python
label_table = qxy.cn_name(adata_cn)
```

The naming rules are applied in order:

- top cell type at least 50%: `<type> hi`;
- top cell type at least 35%: `<type> mix`;
- top two cell types together at least 55%: `<type 1> + <type 2>`; or
- otherwise: `<type 1> + <type 2> mix`.

Duplicate labels receive numeric suffixes. Names are sanitized for safe path
handling on Windows, macOS, and Linux. Shorten long cell-type names before
naming when needed:

```python
label_table = qxy.cn_name(
    adata_cn,
    compaction={"CD8+PD1+LAG3+": "PD1 LAG3 CD8"},
)
```

`qxy.cn_name()` replaces the cluster IDs in `adata_cn.obs["cn"]` with the
descriptive labels and preserves the ID-to-label mapping in
`adata_cn.uns["cn"]["label_map"]`.

## Review and plot neighbourhoods

Inspect the label table and category counts:

```python
label_table
adata_cn.obs["cn"].value_counts(dropna=False)
```

Plot neighbourhoods spatially:

```python
qxy.plot_spatial(
    adata_cn,
    category_col="cn",
    underlay_adata=adata,
    show=True,
)
```

The complete `adata` object can provide a grey underlay when neighbourhoods
were calculated on a filtered copy. Compare neighbourhood frequencies across
samples or groups:

```python
qxy.plot_stacked_bar(adata_cn, category_col="cn")
qxy.plot_stacked_bar(
    adata_cn,
    category_col="cn",
    group_col="group",
)
```

Create a neighbourhood-abundance heatmap across samples:

```python
# CN fractions within each Sample; falls back to Image when Sample is absent.
qxy.plot_cn_heatmap(adata_cn)

# Add an experimental-condition strip when metadata are available.
qxy.plot_cn_heatmap(adata_cn, condition_col="group")

# Also show the sample origin of each neighbourhood.
qxy.plot_cn_heatmap(adata_cn, normalize="both")
```

Neighbourhood classes are descriptive summaries of local composition. Review
their spatial distribution, dominant cell types, sample balance, and stability
across reasonable `k` and `n_cn` choices before biological interpretation.

## Outputs and saving

Neighbourhood analysis adds:

- `adata.obsm["cn_profile"]` with local cell-type fractions;
- `adata.obs["cn"]` with cluster IDs or descriptive labels; and
- `adata.uns["cn"]` with parameters, cell types, and the label map.

`qxy.cn_name()` also writes:

```text
cn/cn_labels.csv
```

Plotting functions write their figures and matrices under `plots/`. Save the
AnnData object after reviewing the analysis:

```python
qxy.save(adata_cn)
```

## Common problems

- **Missing `celltype`:** complete and review cell typing first, or pass the
  correct column through `cell_type_col`.
- **Missing spatial coordinates:** rebuild the AnnData object with
  `qxy.import_cells()`.
- **Missing `cn_profile`:** run `qxy.cn_knn()` before `qxy.cn_kmeans()` or
  `qxy.cn_name()`.
- **Neighbourhoods mix separate images:** keep `image_col="Image"` so KNN
  graphs remain image-specific.
- **Unstable or uninterpretable clusters:** compare reasonable `k`, `n_cn`,
  and random-seed choices and inspect the spatial plots.

## Related documentation

- Prepare reviewed labels with the [cell-typing guide](cell_typing.md).
- Add experimental groups with the [sample-metadata guide](metadata.md).
- See all figure controls in the [plotting guide](plotting.md).
- Review stored arrays and save behavior in
  [AnnData structure and outputs](anndata_and_outputs.md).
