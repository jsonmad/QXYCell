<img src="../assets/qxycell-icon.png" alt="QXYCell logo" width="96">

# QXYCell Overview

A Python bridge from a QuPath project folder to AnnData.

QXYCell turns a QuPath project into an analysis-ready AnnData spatial object. It reads cell level measurement tables and GeoJSON annotation polygons into an `.h5ad` file, then applies marker thresholds and cell type rules as explicit downstream steps for Scanpy, Squidpy, pandas, and custom Python workflows.

**Required preparation:** complete the [QuPath 0.7 multiplex-IF preparation guide](qupath_preparation.md) before running QXYCell. It covers image and pixel-size verification, segmentation, required project-folder assets, filenames, and preflight checks.

## Overview

QXYCell reads one QuPath project folder containing exported measurements, classifier JSON files, GeoJSON annotations, and any project-level threshold tables. Keep all exported inputs inside that folder; QXYCell searches it recursively. It extracts marker measurements, assigns annotation labels from spatial boundaries, and stores the result as an AnnData `.h5ad` object. Marker thresholding and cell typing are separate explicit steps.

The goal is to make QuPath outputs immediately usable with tools like Scanpy, Squidpy, pandas, and custom Python workflows.

## Core Staged Workflow

Stage 1 creates the base measurement checkpoint. Later stages refresh the same active H5AD, `tables/cells_obs.csv`, and `tables/markers_var.csv`. The optional ignore-removal stage deletes cells in annotations drawn around tissue folds, damaged tissue, debris, edge artifacts, staining artifacts, or other excluded regions. Changed ignore polygons therefore require rebuilding stages 1 and 2 first.

    import qxycell as qxy

    qupath_project = "/path/to/qupath_project"

    report = qxy.check(qupath_project)

    # 1. Base AnnData from measurements
    adata = qxy.import_cells(qupath_project)

    # 2. Replace current GeoJSON-derived data
    qxy.add_annotations(adata, pixel_size_um=0.28)

    # Optional 2b. Remove cells in tissue/staining artifact annotations
    qxy.remove_cells(adata, remove_cells="ignore")

    # 3A. Use classifier JSON thresholds only
    qxy.threshold_from_classifiers(adata)
    # Saves/replaces thresholds/classifier_thresholds.tsv

    # Or 3B. Use one named threshold table only
    # qxy.threshold_from_table(adata, "thresholds.tsv")

    # 4. Generate the prompt, then review and save the returned YAML
    qxy.celltype_prompt(adata, context="Describe the tissue and expected populations")

    # 5. Apply the reviewed YAML
    qxy.celltype(adata, "/path/to/celltype_logic.yaml")

    # 6. Optionally save spatial cell-type plots
    qxy.plot_spatial(adata, category_col="celltype", show=False)

### Rerun and replacement rules

| Changed input | Rerun | Replacement behaviour |
|----|----|----|
| Annotation or cell GeoJSON | `qxy.add_annotations(adata)` | Replaces annotation, sample, and cell-polygon columns; invalidates thresholds and later stages. |
| Ignore annotation polygons | Rebuild stages 1 and 2, then run `qxy.remove_cells(adata, remove_cells="ignore")` | Restores the full measurement set before removing cells inside the current tissue or staining artifact regions. |
| Classifier JSON thresholds | `qxy.threshold_from_classifiers(adata)` | Replaces marker positivity and `thresholds/classifier_thresholds.tsv`; ignores supplied threshold tables; invalidates prompt, cell types, and post-analysis. |
| Reviewed threshold table | `qxy.threshold_from_table(adata, table)` | Replaces marker positivity using only the named table; invalidates prompt, cell types, and post-analysis. |
| Prompt context | `qxy.celltype_prompt(adata, context=...)` | Replaces `celltype/current_prompt.txt` but preserves an expert-edited YAML. |
| Cell-type YAML | `qxy.celltype(adata)` | Replaces prior cell type, feature, derived-feature, count, and rule-summary outputs. |

If no output folder is supplied, `qxy.check(...)` writes reports to:

    ../qupath_project_check_YYMMDD_HHMM/

The saved `.h5ad` from `qxy.import_cells(...)` is written by default to:

    ../qupath_project_run_YYMMDD_HHMM/h5ad/qxycell.h5ad

## Inputs

QXYCell requires a manually exported QuPath measurement table. Other inputs are needed only for the corresponding optional feature:

    measurements.csv or measurements.tsv  (required)
    *.json  (threshold option A: QuPath classifier JSON)
    thresholds.tsv or thresholds_YYMMDD-HHMM.tsv  (threshold option B: reviewed table)
    GeoJSON annotation or cell files  (optional spatial annotations and geometry)

The measurement table must contain standard QuPath columns:

    Image
    Object ID
    Centroid X µm
    Centroid Y µm

The known encoding variants `Centroid X ¬µm` and `Centroid Y ¬µm` are also accepted and normalised to the standard coordinate names during import.

QuPath centroid columns are already in micrometres. QuPath GeoJSON is in full-resolution pixel coordinates, so `qxy.add_annotations()` scales it using `pixel_size_um` (default `0.28`). The value must be positive and finite. QXYCell supports square pixels only: verify that QuPath's pixel width and height are equal and pass their single value rather than averaging unequal values.

Thresholding is an explicit choice. `qxy.threshold_from_classifiers()` reads QuPath classifier JSON files only, refuses conflicting definitions, and writes the applied values to the stable `thresholds/classifier_thresholds.tsv` file. `qxy.threshold_from_table()` reads only the named CSV/TSV table and never falls back to classifier JSON. A fresh editable table can be created with `qxy.generate_threshold_table(project_dir, output_dir=output_dir)`. Every successful stage refreshes the active H5AD and exported observation/marker tables.

## Outputs

    qupath_project_check_YYMMDD_HHMM/
      check_report.txt
      check_report.json
      tables/
        measurement_files.csv
        classifier_report.csv
        geojson_report.csv
        validation_messages.csv

    qupath_project_run_YYMMDD_HHMM/
      h5ad/qxycell.h5ad
      celltype/
        current_prompt.txt
        celltype_rules_summary_<logic-name>.tsv
      thresholds/
        classifier_thresholds.tsv
      tables/
        cells_obs.csv
        markers_var.csv
        annotation_assignments.csv
        annotation_conflicts.csv
        thresholding_summary.csv
        celltype_counts.csv
      run.log

The threshold and cell-type files shown above appear after their corresponding stages. `qxy.generate_threshold_table()` adds `thresholds/thresholds_YYMMDD-HHMM.tsv` only when explicitly requested.

## Technical Overview

QXYCell uses separate, rerunnable functions for validation and each core processing stage:

    qxy.check(project_dir)
    adata = qxy.import_cells(project_dir)
    qxy.add_annotations(adata)
    qxy.remove_cells(adata, remove_cells="ignore")  # optional artifact removal
    qxy.threshold_from_classifiers(adata)  # or threshold_from_table(...)
    qxy.celltype_prompt(adata)
    qxy.celltype(adata, "/path/to/celltype_logic.yaml")
    qxy.plot_spatial(adata, category_col="celltype", show=False)

### `check()`

Performs read-only preflight inspection and validation. It discovers measurement tables, threshold definitions, object classifier JSONs, and GeoJSON files; lists every annotation and its planned AnnData assignment; and writes text/JSON reports. It does not generate a threshold table, apply thresholds, apply cell typing, or generate an LLM prompt.

### Staged API

The API separates cell import, annotations, threshold-source selection, prompt generation, and cell typing. Each successful stage checkpoints its owned columns and files so later stages can be rerun explicitly.

## AnnData Layout

`adata.X` contains the raw marker measurement matrix imported by `qxy.import_cells()`. One of the two explicit threshold stages later adds marker positivity columns to `adata.obs`. Re-running a threshold stage replaces its positivity columns and removes downstream cell-type and feature columns because those labels depend on the prior thresholds.

`adata.obs` contains per-cell metadata and computed annotations:

    Image
    Object ID
    Xµm
    Yµm
    quxy_source_file
    quxy_source_row
    quxy_cell_id
    TMA Core              # Optional QuPath measurement-table column
    CoreID                # Categorical — derived only from TMA Core
    <marker>_pos        # added by one explicit threshold stage
    annotation__<label>
    cell_polygon_wkt      # WKT cell boundary polygons when cell GeoJSON is available
    Sample                # Categorical — NaN if unassigned; added by qxy.add_annotations() / qxy.assign_samples()
    celltype              # added by qxy.celltype() / qxy.apply_celltypes()
    cn    # Categorical — added by qxy.cn_kmeans()

Example marker positivity columns:

    CD3_pos
    CD4_pos
    CD8_pos
    PD1_pos
    GZMB_pos

Example annotation columns:

    annotation__Ignore
    annotation__Region
    annotation__Unclassified

`adata.var` contains marker/classifier metadata:

    marker_name
    classifier_name
    source_measurement_column
    threshold
    threshold_source

`adata.obsm["spatial"]` contains cell centroid coordinates from QuPath:

    Xµm
    Yµm

`adata.uns["qxycell"]` contains run metadata, validation status, annotation source-to-column assignments, threshold and cell-type application flags and provenance, and explicit `llm_prompt_generated`/`llm_prompt_path` fields.

## Classifier Handling

In v1, QXYCell applies threshold definitions as a separate step of the form:

    measurement column >= threshold

`qxy.threshold_from_classifiers()` searches recursively for usable single-measurement QuPath classifier JSON files, ignores supplied tables, and writes or replaces `thresholds/classifier_thresholds.tsv` with the values it applied. `qxy.threshold_from_table()` uses one reviewed table, ignores classifier JSON, and does not modify that supplied table. A separate editable table template can be created with `qxy.generate_threshold_table()`. Active threshold definitions are reported in:

    tables/classifier_report.csv

If different classifier JSONs target the same measurement and image scope with different thresholds, QXYCell treats this as an unresolved conflict rather than choosing the last filename. The check audit lists every candidate in `tables/classifier_conflicts.csv`. Generated tables add `classifier_conflict`, `candidate_classifiers`, `candidate_thresholds`, and `candidate_sources` columns, and leave the affected image cells blank.

Import may proceed while the table awaits review, but threshold application is blocked until every image column in each conflict-marked row has a numeric value. Per-image values are applied by exact match to `adata.obs["Image"]`. This supports iterative threshold and cell-type review:

    # Edit one value per image in the generated TSV, then:
    qxy.threshold_from_table(adata, threshold_file)
    qxy.celltype(adata, "celltype_logic.yaml")

    # Edit the same TSV again and repeat both calls.
    # Previous threshold-dependent cell type columns are removed.

## Cell Type Logic Prompt

QXYCell can generate a copy-pasteable prompt for an external LLM to draft a first-pass cell type logic YAML from the markers present in a loaded AnnData object:

    prompt = qxy.celltype_prompt(adata)

By default, the prompt is printed to the active Python session and saved to:

    <run-folder>/celltype/current_prompt.txt

The generated YAML can be saved as:

    <run-folder>/celltype/celltype_logic.yaml

and applied with:

    qxy.celltype(adata, "<run-folder>/celltype/celltype_logic.yaml")

## Cell Neighbourhood Analysis

After cell types are assigned, QXYCell can characterise each cell by the composition of its local spatial neighbourhood using `qxy.cn_knn` and `qxy.cn_kmeans`:

    qxy.cn_knn(adata, k=12)
    qxy.cn_kmeans(adata, n_cn=12)

`cn_knn` builds a k-nearest-neighbour graph per image (to prevent neighbours bleeding across tissue sections), then computes a per-cell composition frequency vector stored as `adata.obsm["cn_profile"]`. `cn_kmeans` clusters those profiles into discrete CN types using MiniBatchKMeans, storing labels `"N0"`…`"N{n-1}"` as a `pd.Categorical` in `adata.obs["cn"]`.

Cell neighbourhoods can be visualised spatially with:

    qxy.plot_spatial(adata, sample_col="ImageID", category_col="cn")

To exclude a cell type entirely, subset a copy before both CN steps. To retain the full tissue as a grey spatial underlay while plotting the filtered CN object:

    adata_cn = adata[adata.obs["celltype"] != "Unknown"].copy()
    qxy.cn_knn(adata_cn, k=12)
    qxy.cn_kmeans(adata_cn, n_cn=12)
    qxy.plot_spatial(adata_cn, underlay_adata=adata,
                     sample_col="Sample", category_col="cn")

## Annotation and TMA Core Handling

QXYCell uses core IDs only when the cell measurement table contains QuPath's exact `TMA Core` column. `qxy.import_cells()` preserves that column and automatically creates categorical `adata.obs["CoreID"]`. Without `TMA Core`, no `CoreID` column is created.

GeoJSON files are also parsed for these feature types:

- **Annotation polygons** (`objectType: annotation`) — regular labels become boolean `annotation__<label>` columns. Labels containing `Sample` collapse into one `adata.obs["Sample"]` column. Ignore annotations remain removable with `qxy.remove_cells()` or `qxy.remove_annotations(adata, text="Ignore")`.
- **Cell polygons** (`objectType: cell`) — matched to cells by QuPath Object ID and stored as WKT strings in `adata.obs["cell_polygon_wkt"]`. These can be plotted with `qxy.plot_cell_boundaries()`.

Cell boundary GeoJSON features are matched to cells by QuPath Object ID. Annotation geometry is applied by spatial containment of cell centroids.

Annotation conflicts are flagged only when the *same annotation column* is triggered more than once for a single cell — i.e. two overlapping polygons of the same label. Overlaps between different annotation labels are not flagged. Conflicts are written to:

    tables/annotation_conflicts.csv

## Current Plotting Defaults

- `plot_spatial()`: automatically uses usable `Sample` labels and otherwise falls back to `Image`; PNG only; cells with missing labels are excluded; a separate full-data underlay can be supplied with `underlay_adata=`.
- `plot_stacked_bar()`: PDF only plus CSV tables.
- `plot_cell_boundaries()`: PNG only.
- `plot_annotation_polygons()`: PNG-only QC, low-resolution cell-density underlay, and boundary-only polygons (`fill=False`).
- `plot_marker_positivity_heatmap()`: fraction positive (0–1), PDF plus CSV.
- `plot_marker_intensity_heatmap()`: category mean intensity followed by per-marker Z-scoring, PDF plus CSV. Default markers are those actually applied by thresholding.
- `plot_cn_heatmap()`: automatically uses usable `Sample` labels and otherwise falls back to `Image`; PDF plus CSV.

## Verification

The public function examples are generated from a small synthetic QuPath project folder. This keeps the documentation reproducible without publishing research data or local filesystem paths.
