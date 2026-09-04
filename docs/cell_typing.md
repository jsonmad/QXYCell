# Cell typing

QXYCell reads marker thresholds from single-measurement classifiers saved in
the QuPath project or from a threshold table and applies them to determine
marker positivity for each cell. QXYCell then generates a prompt for drafting
ordered cell-type rules in YAML. The user—a biology domain expert—reviews the
YAML logic, after which QXYCell applies the reviewed rules to AnnData and
assigns a cell-type label to each cell.

## Prerequisites

Apply exactly one threshold source before cell typing. Use one of the following
routes, not both.

### Apply saved QuPath classifiers

```python
qxy.threshold_from_classifiers(adata)
```

This reads the saved single-measurement classifier JSON files from the QuPath
project, applies their thresholds, and writes the applied values to
`thresholds/classifier_thresholds.tsv` in the active output folder.

### Generate and apply a threshold table

QXYCell can generate a fresh editable table from the measurement columns and
available classifier thresholds in the QuPath project:

```python
threshold_table = qxy.generate_threshold_table(
    "/path/to/qupath_project",
    output_dir="/path/to/output",
)
```

The returned path points to a timestamped TSV in the output folder's
`thresholds/` directory. The table has one threshold column for each image.
Review every row and enter a numeric value for every image before applying it.
Conflict-marked rows may contain blank image values that must be resolved
during this review.

Apply the reviewed table using the returned path:

```python
qxy.threshold_from_table(adata, threshold_table)
```

An existing reviewed CSV or TSV can be supplied directly instead. Table-based
thresholding uses only the named file and does not fall back to the classifier
JSON files.

Thresholding creates `adata.obs["<MARKER>_pos"]` columns. Cell-type rules use
the marker portion of those names exactly, including compartment suffixes such
as `CD3-nuc`, `CD3-mem`, or `CD3-cyto`.

```python
sorted(column for column in adata.obs.columns if column.endswith("_pos"))
```

Review thresholds before writing cell-type rules. A cell-type schema cannot
correct inappropriate positivity calls.

## Generate a first-pass prompt

QXYCell can create a prompt containing the available thresholded markers and
project context:

```python
prompt = qxy.celltype_prompt(
    adata,
    context=(
        "Describe the tissue, disease, experimental groups, and expected "
        "cell populations."
    ),
)
```

The prompt is printed, returned, and saved to:

```text
celltype/current_prompt.txt
```

Copy it into an external LLM to draft `celltype_logic.yaml`. Treat the response
only as a starting point. A biology domain expert must review the marker names,
phenotypes, exclusions, rule order, and fallback behavior before the YAML is
applied.

## Cell-type YAML structure

The file can contain exclusive cell-type `rules`, non-exclusive `features`,
and `derived_features`:

```yaml
rules:
  - name: CD8_T_cell
    positive: [CD3, CD8]
    negative: [CD4]

  - name: CD4_T_cell
    positive: [CD3, CD4]
    negative: [CD8]

  - name: Macrophage
    positive: [CD68]
    any_positive: [CD163, CD206]

features:
  PD1_positive:
    positive: [PD1]

  Proliferating:
    positive: [Ki67]

derived_features:
  PD1_proliferating:
    all_of: [PD1_positive, Proliferating]
```

Rules and features should use marker names present in the generated prompt.
Derived features may also use feature names defined in the same YAML.

## Rule behavior

Cell-type rules are exclusive and evaluated from top to bottom. The first
matching rule assigns `adata.obs["celltype"]`; later matching rules cannot
replace it.

- `positive` requires every listed marker or feature to be positive.
- `negative` requires every listed marker or feature to be negative.
- `any_positive` requires at least one listed marker or feature to be positive.

Place rare and specific phenotypes before broad parent populations. Otherwise,
a broad early rule may capture cells intended for a later subtype. Cells that
match no rule receive `Unknown` by default.

Features are non-exclusive flags and do not assign a cell type. Derived
features can combine markers or previously defined feature columns with:

- `positive` or `all_of` for references that must all be positive;
- `negative` or `none_of` for references that must all be negative; and
- `any_of` for at least one positive reference.

## Apply reviewed logic

Save the reviewed YAML, preferably in the active run folder, and apply it:

```python
summary = qxy.celltype(
    adata,
    "/path/to/output/celltype/celltype_logic.yaml",
)
```

`qxy.apply_celltypes()` is the equivalent long function name. If no YAML path
is supplied, QXYCell uses the newest `.yaml` or `.yml` file in the active
`celltype/` folder and prints the selected path.

Alternative output names can be supplied when required:

```python
summary = qxy.celltype(
    adata,
    "/path/to/celltype_logic.yaml",
    celltype_column="celltype_reviewed",
    unknown_label="Unassigned",
)
```

## Review assignments

Inspect counts and the returned diagnostics before continuing:

```python
adata.obs["celltype"].value_counts(dropna=False)
summary["unknown_count"]
summary["zero_raw_match_rules"]
summary["zero_assigned_rules"]
summary["missing_references"]
summary["multi_rule_match_count"]
```

The rule-summary TSV reports raw matches, assigned cells, cells blocked by
earlier rules, overlap with other rules, and missing references for every
rule. Rules with raw matches but no assignments were completely captured by
earlier rules and usually need their order or definitions reviewed.

Validate the result spatially and against marker patterns:

```python
qxy.plot_spatial(adata, category_col="celltype", show=True)
qxy.plot_marker_positivity_heatmap(
    adata,
    category_col="celltype",
    show=True,
)
qxy.plot_marker_intensity_heatmap(
    adata,
    category_col="celltype",
    show=True,
)
```

These plots are biological sanity checks, not independent validation of the
cell-type schema.

## Outputs and provenance

Cell typing adds or replaces:

- `adata.obs["celltype"]`, or the selected `celltype_column`;
- non-exclusive feature and derived-feature columns from the YAML;
- `adata.uns["qxycell_celltyping"]`; and
- staged provenance in `adata.uns["qxycell"]`.

It also writes:

```text
celltype/celltype_rules_summary_<yaml-name>.tsv
tables/celltype_counts.csv
```

The active H5AD and public observation/marker tables are checkpointed after a
successful application.

## Revise and rerun

Rerunning `qxy.celltype()` replaces columns owned by the previous YAML,
including feature columns no longer present in the new logic. It does not
require reimporting measurements.

If thresholds change, rerun the selected Stage 3 threshold method first, then
regenerate or review the YAML and rerun cell typing. Thresholding removes or
marks stale the dependent prompt, cell-type, and post-analysis outputs.

## Common problems

- **No thresholded markers:** apply Stage 3 before generating the prompt.
- **A rule assigns no cells:** inspect its missing references, raw matches, and
  whether an earlier rule captured the same cells.
- **Too many `Unknown` cells:** check threshold quality, missing phenotypes,
  overly restrictive exclusions, and unsupported marker combinations.
- **Unexpected broad population:** move specific rules earlier or add
  biologically justified negative markers to the broad rule.
- **Wrong YAML selected automatically:** pass the intended YAML path
  explicitly.

## Related documentation

- Prepare and revise thresholds with
  [QuPath inputs, annotations, and thresholds](qupath_inputs.md).
- Follow checkpoint and rerun order in
  [Running the staged workflow](running_qxycell.md).
- Create validation figures with the [plotting guide](plotting.md).
- Continue to [cellular neighbourhood analysis](cellular_neighbourhoods.md).
