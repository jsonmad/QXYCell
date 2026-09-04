# Running the staged QXYCell workflow

QXYCell is designed as a reviewable sequence of stages. Each successful stage
updates the active AnnData object and its saved checkpoint, so thresholds,
annotations, prompts, and cell-type logic can be revised independently.

Before starting, prepare the files described in the
[QuPath preparation guide](qupath_preparation.md) and
[input guide](qupath_inputs.md).

## Preflight check

Validate the project folder without running the analysis:

```python
import qxycell as qxy

report = qxy.check("/path/to/qupath_project")
```

`qxy.check()` writes a timestamped sibling folder:

```text
qupath_project_check_YYMMDD_HHMM/
```

The report describes discovered measurement tables, GeoJSON files,
annotations, classifiers, threshold files, conflicts, and expected AnnData
columns. It separates definitions from actions: no thresholds, cell typing, or
LLM prompt generation are performed.

## Core staged workflow

```python
import qxycell as qxy

# Stage 1: import cell measurements and create AnnData.
adata = qxy.import_cells("/path/to/qupath_project")

# Stage 2: add or refresh annotations and cell polygons.
qxy.add_annotations(adata, pixel_size_um=0.28)

# Optional Stage 2b: remove cells in matching artifact annotations.
qxy.remove_cells(adata, remove_cells="ignore")

# Stage 3: choose exactly one threshold source.
qxy.threshold_from_classifiers(adata)
# qxy.threshold_from_table(adata, "/path/to/thresholds.tsv")

# Stage 4: generate the prompt used to draft celltype_logic.yaml.
qxy.celltype_prompt(
    adata,
    context="Describe the tissue and expected populations",
)

# Pause for expert review, then save the returned YAML.
# Stage 5: assign cell types from the reviewed YAML.
qxy.celltype(adata, "/path/to/celltype_logic.yaml")

# Stage 6: plot assigned cell types.
qxy.plot_spatial(adata, category_col="celltype", show=True)
```

Run this sample in one persistent `ipython` or `python` session so `adata`
remains available between stages.

## Checkpoints and active output folder

Stage 1 creates the AnnData object and saves the initial H5AD measurement
checkpoint. Later stages update that same H5AD and refresh:

- `tables/cells_obs.csv`
- `tables/markers_var.csv`

`qxy.import_cells()` creates the timestamped run folder used by downstream
stages:

```text
qupath_project_run_YYMMDD_HHMM/
```

The active path is stored in `adata.uns["qxycell"]["output_dir"]`.
Annotation import, thresholding, metadata, cell typing, plotting, and
`qxy.save()` reuse it when `output_dir` is omitted.

Typical output structure:

```text
qupath_project_run_YYMMDD_HHMM/
  h5ad/
    qxycell.h5ad
  tables/
    cells_obs.csv
    markers_var.csv
  thresholds/
  celltype/
  dataset_summary/
  plots/
  run.log
```

## What to rerun after changes

| Changed input | Rerun | What is replaced |
|---|---|---|
| Annotation or cell GeoJSON | `qxy.add_annotations(adata)` | Annotation, sample, and cell-polygon columns; downstream stages become stale |
| Ignore annotation polygons | Stages 1 and 2, then `qxy.remove_cells(adata, remove_cells="ignore")` | All cells are rebuilt before cells in matching regions are removed |
| QuPath classifier JSON thresholds | `qxy.threshold_from_classifiers(adata)` | Marker `_pos` columns and `thresholds/classifier_thresholds.tsv`; prompt, cell types, and post-analysis become stale |
| Reviewed threshold table | `qxy.threshold_from_table(adata, table)` | Marker `_pos` columns; prompt, cell types, and post-analysis become stale |
| Cell-type YAML | `qxy.celltype(adata)` | `celltype`, feature, derived-feature, count, and rule-summary outputs |

Changed annotation or cell GeoJSON does not normally require reimporting the
measurement table. Ignore-region changes are the exception because cells may
already have been deleted. Reimport measurements first to restore the full
cell set, then refresh annotations and repeat removal.

The two threshold functions are deliberately separate. Classifier mode ignores
tables; table mode uses only the named table. Each classifier-mode run replaces
`thresholds/classifier_thresholds.tsv` with the values actually applied.
To refine those values in table mode, first copy or rename that file so a later
classifier-mode run cannot overwrite the manual edits. See
[Refine Stage 3A thresholds with Stage 3B](qupath_inputs.md#refine-stage-3a-thresholds-with-stage-3b).

## Annotation and measurement outputs

After Stage 2, the main AnnData locations are:

| Location | Contents |
|---|---|
| `adata.obs` | Per-cell `Image`, `Object ID`, `Xµm`, `Yµm`, optional `TMA Core` and `CoreID`, `annotation__<label>` columns, optional `Sample`, and optional `cell_polygon_wkt` |
| `adata.X` | Marker intensity matrix, cells × markers |
| `adata.var` | Marker names and metadata |
| `adata.obsm["spatial"]` | Cell centroid coordinates in micrometres |
| `adata.uns["qxycell"]` | Run metadata, output paths, stage status, and provenance |

See [AnnData and outputs](anndata_and_outputs.md) for the complete data model.

## Explicit pixel size

Supply the verified square-pixel size when it differs from 0.28 µm:

```python
qxy.add_annotations(adata, pixel_size_um=0.325)
```

The value scales annotation and cell GeoJSON coordinates from full-resolution
pixels into micrometres. QuPath centroid measurements are already in
micrometres.

## Optional single-call workflow

The staged API is recommended for reviewable, rerunnable analyses.
`qxy.workflow()` remains available for a single convenience call:

```python
adata = qxy.workflow(
    "/path/to/qupath_project",
    threshold_file="/path/to/reviewed_thresholds.tsv",
    apply_thresholds=True,
    sample_metadata="sample_metadata.tsv",
    sample_col="Image",
    celltype_logic="celltype_logic.yaml",
)
```

## Continue after the staged workflow

- Add experimental fields with the [sample-metadata guide](metadata.md).
- Create and review assignments with the [cell-typing guide](cell_typing.md).
- Analyse local composition with the
  [cellular-neighbourhood guide](cellular_neighbourhoods.md).
- Create figures with the [plotting guide](plotting.md).
- Save, reload, and inspect the object with
  [AnnData and outputs](anndata_and_outputs.md).
