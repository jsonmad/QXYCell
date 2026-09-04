# QXYCell technical workflow

A staged QuPath-to-AnnData workflow with explicit checkpoints, a mutually exclusive threshold choice, and targeted rerun dependencies.

**Legend:** Optional = dashed treatment in the original diagram; Decision = choose exactly one path; Domain expert review = human approval gate.

## Forward workflow

Stages 1–5 update the active H5AD and refresh both `tables/cells_obs.csv` and `tables/markers_var.csv`. Optional Stage 2b refreshes the filtered H5AD and `cells_obs.csv`; preflight and Stage 6 create reports or plots without changing the active checkpoint.

1.  **Prep — QuPath**

    ### Prepare QuPath project

    Export the measurements and spatial assets that define the run inputs.

    Prepare the measurement table (`measurements.tsv`), annotation GeoJSON (`slide01.geojson`), cell segmentation GeoJSON (`slide01-cells.geojson`), and threshold assets (`thresholds.tsv` or `classifiers/object_classifiers/*.json`) before import.

    **Optional preflight:** Run `qxy.check(project_dir)` to validate the QuPath project and write a timestamped report before creating an analysis checkpoint.

2.  **Stage 1 — Python · QXYCell**

    ### Stage 1 — `qxy.import_cells(project_dir)`

    Create the base AnnData object from the QuPath measurement table.

    Populates `adata.X` with marker intensities, `adata.obs` with cell identifiers and centroid coordinates (plus `TMA Core` and `CoreID` when available), `adata.var` with marker metadata, `adata.obsm["spatial"]` with x–y coordinates, and `adata.uns["qxycell"]` with run provenance and output paths.

3.  **Stage 2 — Python · QXYCell**

    ### Stage 2 — `qxy.add_annotations(adata)`

    Add or refresh GeoJSON-derived annotations, sample assignments, and cell polygons.

    Annotations with `Sample` in their label create `adata.obs["Sample"]`; other labels create Boolean `adata.obs["annotation__<label>"]` columns; cell segmentation GeoJSON adds `adata.obs["cell_polygon_wkt"]`. GeoJSON coordinates are converted from pixels to micrometres—confirm `pixel_size_um` matches the image calibration (default: `0.28` µm/pixel).

4.  **Optional Stage 2b — Python · QXYCell**

    ### Optional Stage 2b — `qxy.remove_cells(adata, remove_cells="<label>")`

    Remove cells where a matching `adata.obs["annotation__<label>"]` column is `True`.

    Filters the AnnData object in place by removing cells located within the indicated annotation polygons. `<label>` is matched case-insensitively against annotation column names; the default is `"ignore"`.

5.  **Stage 3 — Python · QXYCell — choose one route**

    ### Stage 3 — Apply marker thresholds

    Choose exactly one route:

    #### 3A · Apply classifier JSON directly

        qxy.threshold_from_classifiers(
            adata
        )

    Applies the QuPath classifier thresholds and saves the applied values to `thresholds/classifier_thresholds.tsv`.

    #### 3B · Generate and review a threshold table

        threshold_file = qxy.generate_threshold_table(
            project_dir
        )
        # Review or edit the generated table
        qxy.threshold_from_table(
            adata, threshold_file
        )

    Builds a timestamped threshold table from the classifier JSON files, then applies only the reviewed table.

    3A and 3B are alternatives; run one, not both.

6.  **Stage 4 — Python · QXYCell + LLM**

    ### Stage 4 — `qxy.celltype_prompt(adata, context="<biological context>")`

    Create a project-specific LLM prompt for drafting cell-type logic.

    Builds the prompt from the available `adata.obs["<marker>_pos"]` columns and the supplied biological context. It prints and returns the prompt and saves it as `celltype/current_prompt.txt`. Copy the prompt into an LLM, use its response as the first draft of `celltype_logic.yaml`, then review and correct the rules before saving the file.

7.  **Review gate — Domain expert · YAML**

    ### Domain expert review — review and update `celltype_logic.yaml`

    A domain expert should verify and edit the LLM-generated logic before cell typing.

    Cell-type `rules` are evaluated from top to bottom, and the first matching rule assigns `adata.obs["celltype"]`; unmatched cells remain `Unknown`. Each rule can require positive and negative marker calls, so specific populations should appear before broader populations. Optional `features` and `derived_features` create additional non-exclusive columns. Update the YAML after review and save the corrected version as `celltype/celltype_logic.yaml` inside the active QXYCell run folder.

8.  **Stage 5 — Python · QXYCell**

    ### Stage 5 — `qxy.celltype(adata)`

    Apply the domain-expert-reviewed YAML to assign cell types.

    Writes assignments to `adata.obs["celltype"]`, creates any configured feature and derived-feature columns, saves `tables/celltype_counts.csv` and a rule-summary TSV, and records diagnostics in `adata.uns["qxycell_celltyping"]`.

9.  **Optional Stage 6 — Python · QXYCell**

    ### Optional Stage 6 — Plot assigned cell types spatially

        # One plot per Sample annotation
        qxy.plot_spatial(
            adata,
            sample_col="Sample",
            show=False,
        )

        # One plot per QuPath image
        qxy.plot_spatial(
            adata,
            sample_col="Image",
            show=False,
        )

    Create spatial plots of the assigned cell types, grouped by `Sample` or `Image`.

    If `sample_col` is omitted, QXYCell uses `Sample` when usable sample labels exist and otherwise falls back to `Image`. `show=False` is recommended for scripted runs; PNG output is enabled by default.

10. **Downstream — Python analysis**

    ### Downstream analysis — Scanpy, Squidpy, scimap, pandas, or custom Python

    Load the completed AnnData object from the default sibling run folder:

        qupath_project_run_YYMMDD_HHMM/h5ad/qxycell.h5ad

        adata = qxy.load("path/to/qxycell.h5ad")

    **Important contents**

    - `adata.X`: cells × markers intensity matrix.
    - `adata.var`: marker index with default columns `marker_name`, `classifier_name`, `source_measurement_column`, `threshold`, and `threshold_source`.
    - `adata.obs`: default import columns `Image`, `Object ID`, `Xµm`, `Yµm`, `quxy_source_file`, `quxy_source_row`, and `quxy_cell_id`.
    - Optional or project-derived `adata.obs` columns: `TMA Core`, `CoreID`, `Sample`, `annotation__<label>`, `cell_polygon_wkt`, `<marker>_pos`, `celltype`, and YAML-defined feature columns.
    - `adata.obsm["spatial"]`: cell x–y coordinates in micrometres.
    - `adata.uns["qxycell"]`: run paths, provenance, and stage status.
    - Additional summaries may include `qxycell_annotation_labels`, `qxycell_sample_annotations`, `qxycell_core_ids_from_measurements`, `qxycell_thresholding`, and `qxycell_celltyping`.

    The saved H5AD contains the completed QXYCell pipeline state and can be loaded directly into Scanpy, Squidpy, scimap, pandas, or custom Python for clustering, spatial analysis, cell–cell interaction testing, and cellular-neighbourhood analysis.

## Rerun dependencies

- **Annotation/cell GeoJSON changed** → Stage 2, optional Stage 2b, then Stages 3–5

- **Ignore polygons changed** → rebuild Stages 1–2, then Stage 2b and Stages 3–5

- **Thresholds changed** → selected Stage 3 path, then Stages 4–5

- **Prompt context changed** → Stage 4; preserve expert-edited YAML

- **Cell-type YAML changed** → Stage 5 and optional Stage 6

Supporting documentation: [QuPath preparation](qupath_preparation.md) · [QXYCell overview](QXYCell_overview.md) · [Function reference](QXYCell_function_reference.md).
