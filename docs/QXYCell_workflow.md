# QXYCell technical workflow

A staged QuPath-to-AnnData workflow with explicit checkpoints, a final threshold source, and targeted rerun dependencies.

**Legend:** Optional = dashed treatment in the original diagram; Decision = choose the final threshold source; Domain expert review = human approval gate.

## Forward workflow

Stages 1–5 update the active H5AD and refresh both `tables/cells_obs.csv` and `tables/markers_var.csv`. Optional Stage 2b refreshes the filtered H5AD and `cells_obs.csv`; Stage 6 creates plots without changing the active checkpoint.

## Prep — QuPath

### Prepare QuPath project

Export `measurements.tsv`, annotation GeoJSON (`slide01.geojson`), optional cell segmentation GeoJSON (`slide01-cells.geojson`), and either classifier JSON files (`classifiers/object_classifiers/CD3.json`) or a reviewed threshold table (`thresholds.tsv`) before import.

---

## Stage 1 — Python · QXYCell

### `qxy.import_cells(project_dir)`

Create the base AnnData object from the QuPath measurement table. It stores marker intensities in `adata.X`, cell identifiers and centroids in `adata.obs`, spatial coordinates in `adata.obsm["spatial"]`, and run provenance in `adata.uns["qxycell"]`.

---

## Stage 2 — Python · QXYCell

### `qxy.add_annotations(adata)`

Add or refresh GeoJSON-derived annotations, sample assignments, and optional cell polygons. Verify `pixel_size_um` against image calibration.

---

## Optional Stage 2b — Python · QXYCell

### `qxy.remove_cells(adata, remove_cells="<label>")`

Remove cells inside matching annotation polygons; the default label is `"ignore"`.

---

## Stage 3 — Python · QXYCell — 3A OR 3B

### Apply marker thresholds

Use 3A to apply classifier JSON values directly, or use 3B to apply a reviewed threshold table.

#### 3A · Apply classifier JSON directly

    qxy.threshold_from_classifiers(adata)

This writes per-image values to `thresholds/classifier_thresholds.tsv`. Copy or rename it before manually adjusting images whose classifier-derived threshold needs changing, then use that reviewed copy through 3B.

#### 3B · Generate or refine a threshold table

    threshold_file = qxy.generate_threshold_table(project_dir)
    # Review or edit the generated table
    qxy.threshold_from_table(adata, threshold_file)

Build a timestamped table, or use the reviewed copy from 3A. Review and save every marker and image value, then run `qxy.threshold_from_table(adata, threshold_file)` using that reviewed table.

---

## Stage 4 — Python · QXYCell + LLM

### `qxy.celltype_prompt(adata, context="<biological context>")`

Creates a project-specific LLM prompt for drafting cell-type logic. Review and correct the returned YAML before use.

---

## Review gate — Domain expert · YAML

### Review and update `celltype_logic.yaml`

A domain expert should verify marker names, phenotypes, exclusions, rule order, and fallback behaviour before cell typing.

---

## Stage 5 — Python · QXYCell

### `qxy.celltype(adata)`

Apply the domain-expert-reviewed YAML to assign cell types.

---

## Optional Stage 6 — Python · QXYCell

### Plot assigned cell types spatially

    qxy.plot_spatial(adata, sample_col="Sample", show=False)

Create spatial plots grouped by `Sample` or `Image`.

---

## Downstream — Python analysis

### Scanpy, Squidpy, scimap, pandas, or custom Python

Load the completed AnnData object from `qupath_project_run_YYMMDD_HHMM/h5ad/qxycell.h5ad` with `qxy.load()`.

## Rerun dependencies

- **Annotation/cell GeoJSON changed** → Stage 2, optional Stage 2b, then Stages 3–5

- **Ignore polygons changed** → rebuild Stages 1–2, then Stage 2b and Stages 3–5

- **Thresholds changed** → selected Stage 3 path, then Stages 4–5

- **Prompt context changed** → Stage 4; preserve expert-edited YAML

- **Cell-type YAML changed** → Stage 5 and optional Stage 6

Supporting documentation: [QuPath preparation](qupath_preparation.md) · [QXYCell overview](QXYCell_overview.md) · [Function reference](QXYCell_function_reference.md).
