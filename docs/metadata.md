# Sample metadata

QXYCell can attach sample-level experimental metadata to every matching cell in
an AnnData object. Typical fields include treatment group, patient or animal
identifier, time point, tissue site, batch, and outcome. Metadata import does
not change marker measurements, annotations, thresholds, or cell types.

## Shorten image names

Use `qxy.simple_image_names()` to retain the text after the last underscore
and remove a trailing `.ome.tiff`:

```python
mapping = qxy.simple_image_names(adata)
adata.obs[["Image_original", "Image"]].drop_duplicates()
```

For example, `scan_Buckingham_21P10459A1.ome.tiff` becomes `21P10459A1`.
The function updates `Image` in place and copies the original column to
`Image_original` on its first call. Subsequent calls use that preserved column.
Distinct images with the same short name receive `_2`, `_3`, and so on, in
first-appearance order; cells belonging to the same image share one name.
Missing values remain missing. The returned dictionary maps original names to
short names. Only these two observation columns are changed.

Annotation assignment and default image-specific threshold matching use
`Image_original` when present, so existing QuPath exports still match after
renaming. Annotation plots also match geometry and cell underlays using the
original names while displaying the short names. Keep `Image_original` intact.

By default, changes stay in memory. To overwrite the loaded H5AD with the
updated AnnData object:

```python
adata = qxy.load("my_data.h5ad")
mapping = qxy.simple_image_names(adata, save=True)
```

`save=True` writes the entire current AnnData object, including other unsaved
changes. CSV files are not updated. The save path must be known and exist;
use `qxy.load()` (or `qxy.load_latest()`) to record the file actually opened.
For an object loaded with another library, first use
`qxy.save(adata, "my_data.h5ad")` to establish the intended save path.

When matching metadata that still contains the full image names, use
`sample_col="Image_original"` with `qxy.add_metadata()`.

## Prepare the metadata table

Use a CSV or TSV with one row per sample and one column containing the sample
identifier. The remaining columns contain values to add to `adata.obs`.

```text
sample    group      timepoint    batch
Sample_A  control    baseline     1
Sample_B  treated    day_7        1
Sample_C  treated    day_14       2
```

Each identifier must be unique. Matching is case-sensitive after surrounding
whitespace is removed. Choose an AnnData identifier column whose values refer
to the same units as the metadata rows:

- `Image` for one metadata row per QuPath image.
- `Sample` for regions defined by annotations containing `Sample` in their
  names.
- `CoreID` for one metadata row per imported TMA core.

Inspect the available identifiers before importing metadata:

```python
adata.obs["Sample"].value_counts(dropna=False)
```

## Add metadata

Specify both matching columns explicitly when their names differ:

```python
summary = qxy.add_metadata(
    adata,
    "sample_metadata.tsv",
    sample_col="Sample",
    metadata_sample_col="sample",
)
```

`sample_col` identifies the key in `adata.obs`, while
`metadata_sample_col` identifies the corresponding key in the CSV or TSV.
When `metadata_sample_col` is omitted, QXYCell uses a table column matching
`sample_col`, or otherwise the first table column. Explicit names are
recommended because they make the match unambiguous.

The imported values are broadcast to every cell with the matching sample
identifier:

```python
adata.obs[["Sample", "group", "timepoint", "batch"]].head()
```

## Select or rename imported columns

Import only selected fields:

```python
qxy.add_metadata(
    adata,
    "sample_metadata.tsv",
    sample_col="Sample",
    metadata_sample_col="sample",
    columns=["group", "timepoint"],
)
```

Add a prefix when the source names could conflict with existing AnnData
columns:

```python
qxy.add_metadata(
    adata,
    "sample_metadata.tsv",
    sample_col="Sample",
    metadata_sample_col="sample",
    prefix="clinical_",
)
```

QXYCell refuses to replace an existing `adata.obs` column unless
`overwrite=True` is supplied. Use overwriting only after confirming that the
new table is the intended source:

```python
qxy.add_metadata(
    adata,
    "sample_metadata.tsv",
    sample_col="Sample",
    metadata_sample_col="sample",
    columns=["group"],
    overwrite=True,
)
```

A pandas DataFrame can be supplied instead of a file path.

## Review matching results

`qxy.add_metadata()` returns a summary containing:

- the AnnData and metadata matching columns;
- the number of identifiers in each source;
- the number of matched identifiers;
- AnnData identifiers missing from the metadata table;
- unused metadata rows; and
- columns added to `adata.obs`.

```python
summary["n_matched_samples"]
summary["missing_in_metadata"]
summary["unused_metadata"]
```

Unmatched cells receive missing values in the imported columns. Duplicate keys
in the metadata table raise an error because QXYCell cannot determine which
row should be used.

## Outputs and saving

The matching summary is stored in
`adata.uns["qxycell_sample_metadata"]`. QXYCell also writes:

```text
tables/sample_metadata_summary.tsv
tables/sample_metadata_applied.tsv
```

The first table records matching results. The second preserves the metadata
table that was applied. Save the AnnData checkpoint after reviewing the match:

```python
qxy.save(adata)
```

## Common problems

- **No matches:** compare the exact values in the two identifier columns and
  confirm that they represent the same unit, such as image versus sample.
- **Some samples are missing:** inspect `missing_in_metadata` and
  `unused_metadata` for spelling, capitalization, or naming differences.
- **Duplicate-key error:** reduce the metadata table to one row per identifier
  before importing it.
- **Column already exists:** use a prefix, choose a different source column,
  or deliberately set `overwrite=True`.

## Related documentation

- See how `Image`, `Sample`, and `CoreID` are created in
  [QuPath inputs, annotations, and thresholds](qupath_inputs.md).
- Review stored fields and checkpoints in
  [AnnData structure and outputs](anndata_and_outputs.md).
- Use imported experimental groups in the [plotting guide](plotting.md).
