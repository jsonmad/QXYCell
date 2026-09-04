<img src="assets/qxycell-icon.png" alt="QXYCell logo" width="96">

# QXYCell Function and Command Reference

Python and CLI reference for converting one QuPath project folder into an AnnData workflow.

**Standard import:** `import qxycell as qxy`

**Prepare first:** follow the [QuPath 0.7 multiplex-IF preparation guide](qupath_preparation.md) to verify calibration and create the required project-folder assets.

**Recommended staged workflow:** `import_cells → add_annotations → optional remove_cells → threshold_from_classifiers/table → celltype_prompt → celltype → optional plot_spatial`. Stages replace their outputs in the same active H5AD and exported tables. If ignore polygons change after cells have been removed, rebuild stages 1 and 2 before running the optional removal stage again.

## Core QXYCell Functions

| Function | Purpose | Typical command | Key output |
|----|----|----|----|
| `qxy.check` | Inspect and validate a QuPath project folder. The report lists all annotation names and planned AnnData assignments, identifies the active threshold-definition source, and explicitly records that check does not apply thresholds, apply cell typing, or generate an LLM prompt. | `report = qxy.check(qupath_project)` | `CheckReport` plus files in `<project-name>_check_YYMMDD_HHMM/`. |
| `qxy.generate_threshold_table` | Create a fresh timestamped threshold table without modifying existing files. Competing thresholds for the same measurement are preserved as candidates while affected per-image cells are left blank. | `path = qxy.generate_threshold_table(qupath_project)` | `thresholds_YYMMDD-HHMM.tsv` under the output `thresholds/` folder. |
| `qxy.import_cells` | Stage 1: create the base AnnData from measurement CSV/TSV files only. Marker intensities are stored in `adata.X`; cell identity and coordinates are stored in `adata.obs` and `adata.obsm["spatial"]`. | `adata = qxy.import_cells(qupath_project)` | Base H5AD plus synchronized observation and marker tables. |
| `qxy.add_annotations` | Stage 2: discover current GeoJSON files and replace GeoJSON-derived annotation, sample, and cell-polygon columns. Rerunning invalidates downstream stages. | `summary = qxy.add_annotations(adata, pixel_size_um=0.28)` | Refreshed annotation columns, H5AD, and `cells_obs.csv`. |
| `qxy.threshold_from_classifiers` | Stage 3A: replace marker positivity using QuPath classifier JSON thresholds only. Threshold tables are ignored, conflicting classifiers cause an error, and the applied values are saved to a stable threshold table. | `summary = qxy.threshold_from_classifiers(adata)` | Refreshed `<marker>_pos` columns, H5AD, tables, stage record, and `thresholds/classifier_thresholds.tsv`. |
| `qxy.threshold_from_table` | Stage 3B: replace marker positivity using one named threshold table only. Classifier JSON thresholds are not used as a fallback. | `summary = qxy.threshold_from_table(adata, "thresholds.tsv")` | Refreshed `<marker>_pos` columns, H5AD, tables, and stage record. |
| `qxy.threshold` | Legacy automatic-source threshold entry point retained for compatibility. New workflows should select one of the two explicit threshold functions. | `summary = qxy.threshold(adata, qupath_project)` | `adata.obs["<marker>_pos"]` columns and `adata.uns["qxycell_thresholding"]`. |
| `qxy.apply_thresholds` | Long-name alias for `qxy.threshold`. | `summary = qxy.apply_thresholds(adata, threshold_file="thresholds.tsv")` | The same positivity columns and thresholding summary as `qxy.threshold`. |
| `qxy.celltype_prompt` | Stage 4: create a copy-paste LLM prompt from thresholded markers. Rerunning replaces `celltype/current_prompt.txt` but never overwrites an expert-edited YAML. | `prompt = qxy.celltype_prompt(adata)` | Prompt text and a saved prompt file in the run `celltype/` folder. |
| `qxy.celltype` | Stage 5: replace cell type, feature, and derived-feature columns using ordered first-match YAML logic. | `summary = qxy.celltype(adata)` | Adds `adata.obs["celltype"]` and saves a rule summary TSV. |
| `qxy.apply_celltypes` | Long-name alias for `qxy.celltype`. | `summary = qxy.apply_celltypes(adata, "celltype_logic.yaml")` | Adds `adata.obs["celltype"]` and saves a rule summary TSV. |
| `qxy.remove_cells` | Remove cells in annotation columns identifying tissue folds, damaged tissue, debris, edge artifacts, staining artifacts, or other excluded regions. | `qxy.remove_cells(adata, remove_cells="ignore")` | Filters `adata` in place by default. |
| `qxy.remove_annotations` | Remove cells inside annotation columns matching a custom text string. | `qxy.remove_annotations(adata, text="artifact")` | Filters `adata` in place by default. |
| `qxy.assign_samples` | Assign one sample label per cell from annotation columns containing `sample`. The result is a `pd.Categorical`: cells in exactly one sample annotation receive that label; cells in multiple annotations receive `"Ambiguous"`; cells outside all sample annotations receive `NaN` and are excluded from `value_counts()` and `groupby()`. | `sample_summary = qxy.assign_samples(adata)` | Adds `adata.obs["Sample"]` and warns on overlapping sample annotations. |
| `qxy.assign_annotations` | Collapse selected boolean annotation columns into one categorical observation column. | `summary = qxy.assign_annotations(adata, ["Tumour", "Stroma"], target_col="Region")` | Adds `adata.obs["Region"]` and reports overlaps and unassigned cells. |
| `qxy.assign_core_ids_from_measurements` | Create categorical CoreID from QuPath's measurement-table `TMA Core` column. | `coreid = qxy.assign_core_ids_from_measurements(adata)` | Adds categorical `adata.obs["CoreID"]` from `TMA Core`. Parent and annotations are not accepted as sources. |
| `qxy.add_metadata` | Add sample-level metadata from CSV, TSV, or DataFrame to `adata.obs`. | `qxy.add_metadata(adata, "metadata.tsv", sample_col="ImageID")` | Adds metadata columns to every cell by sample. |
| `qxy.dataset_summary` | Create descriptive dataset-summary tables and an HTML report. | `summary = qxy.dataset_summary(adata, sample_col="Sample")` | Dataset-summary TSV tables and `dataset_summary.html`. |
| `qxy.plot_stacked_bar` | Plot stacked cell type frequencies per sample or group. | `qxy.plot_stacked_bar(adata, sample_col="ImageID")` | PDF plus frequency CSV tables by default; optional PNG. |
| `qxy.plot_spatial` | Plot spatial category maps with shared scale, centered samples, legend, scale bar, optional separate AnnData underlay, and missing-sample exclusion by default. | `qxy.plot_spatial(adata, sample_col="ImageID")` | PNG by default; optional PDF and combined plot. |
| `qxy.plot_annotation_polygons` | Reload original QuPath annotation polygons per image with a low-resolution cell-density underlay. | `qxy.plot_annotation_polygons(adata)` | PNG-only QC plots; boundary-only polygons by default. |
| `qxy.load_cell_polygons` | Load cell boundary polygons from matching QuPath GeoJSON files into an existing AnnData object. | `n = qxy.load_cell_polygons(adata, qupath_project)` | Adds WKT geometry in `adata.obs["cell_polygon_wkt"]`. |
| `qxy.plot_marker_positivity_heatmap` / `qxy.plot_marker_intensity_heatmap` | Two explicit quantities: fraction positive from `_pos` columns, or mean intensity from `adata.X` Z-scored per marker. | `qxy.plot_marker_positivity_heatmap(adata)` or `qxy.plot_marker_intensity_heatmap(adata)` | One PDF plus CSV by default; optional PNG, SVG, or TIFF. |
| `qxy.plot_marker_heatmap` | Compatibility entry point for older code. New code should call the explicit positivity or intensity heatmap function. | `qxy.plot_marker_heatmap(adata, values="positivity")` | The corresponding marker heatmap and CSV matrix. |
| `qxy.plot_cn_heatmap` | Heatmap of CN abundance across samples, normalised by sample or by CN. | `qxy.plot_cn_heatmap(adata)` | PDF plus CSV by default; optional PNG, SVG, or TIFF. |
| `qxy.save` | Save the current AnnData object to the run H5AD path. | `qxy.save(adata)` | Saved `qxycell.h5ad`. |
| `qxy.load` | Load a saved QXYCell H5AD from an H5AD path or output folder. | `adata = qxy.load("<project-name>_run_YYMMDD_HHMM")` | `AnnData`. |
| `qxy.load_latest` | Load the latest timestamped QXYCell run in a folder. | `adata = qxy.load_latest()` | `AnnData`. |
| `qxy.workflow` | Convenience wrapper for run, metadata, ignore removal, celltyping, dataset summary, plots, and save. | `adata = qxy.workflow(qupath_project)` | Completed AnnData workflow. |

## Helper Functions

| Function | Purpose | Command |
|----|----|----|
| `qxy.find_latest_celltype_yaml` | Find the newest YAML file in a celltype folder. | `yaml_path = qxy.find_latest_celltype_yaml()` |
| `qxy.load_celltype_logic` | Load a celltype YAML file as a Python dictionary. | `logic = qxy.load_celltype_logic("celltype_logic.yaml")` |
| `qxy.CheckReport` | Dataclass returned by `qxy.check()`. | `report.n_errors`, `report.n_warnings`, `report.messages` |

## Common Python Workflow

    import qxycell as qxy

    qupath_project = "/path/to/qupath_project"

    adata = qxy.import_cells(qupath_project)
    qxy.add_annotations(adata)
    qxy.remove_cells(adata, remove_cells="ignore")  # optional artifact removal
    qxy.threshold_from_classifiers(adata)  # or threshold_from_table(...)
    # Classifier mode saves/replaces thresholds/classifier_thresholds.tsv

    prompt = qxy.celltype_prompt(
        adata,
        context="Add biological/project context here."
    )

    # Save and review the YAML returned by an LLM, then:
    summary = qxy.celltype(adata, "/path/to/celltype_logic.yaml")
    qxy.plot_spatial(adata, category_col="celltype", show=False)

    qxy.add_metadata(
        adata,
        "sample_metadata.tsv",
        sample_col="ImageID",
    )

    qxy.dataset_summary(adata, sample_col="Sample")

    qxy.plot_spatial(
        adata,
        sample_col="ImageID",
        combined=True,
        max_cols=2,
    )

    qxy.plot_stacked_bar(
        adata,
        sample_col="ImageID",
        group_col="group",
    )

    qxy.save(adata)

## Metadata Commands

| Goal | Command | Notes |
|----|----|----|
| Add metadata matched on QuPath image names | `qxy.add_metadata(adata, "metadata.tsv")` | Default `sample_col="Image"`. |
| Add metadata matched on short sample names | `qxy.add_metadata(adata, "metadata.tsv", sample_col="ImageID")` | `ImageID` must already exist in `adata.obs`. |
| Metadata file uses a different key column | `qxy.add_metadata(adata, "metadata.tsv", sample_col="ImageID", metadata_sample_col="sample_id")` | Use when metadata has `sample_id` but AnnData has `ImageID`. |
| Add only selected metadata columns | `qxy.add_metadata(adata, "metadata.tsv", columns=["group", "mouse_id"])` | Prevents importing unused metadata columns. |
| Overwrite existing metadata columns | `qxy.add_metadata(adata, "metadata.tsv", overwrite=True)` | Use deliberately. |

## Classifier Conflict Review

QXYCell defines a conflict as two or more distinct numeric thresholds for the same measurement column and the same image scope. Different explicitly image-scoped thresholds are valid. Global classifier JSON conflicts are never resolved by filename order.

| Stage | Behaviour |
|----|----|
| `qxy.check()` | Reports each measurement, classifier name, threshold, and source. Writes `tables/classifier_conflicts.csv`. |
| Table generation | Adds candidate-provenance columns, sets `classifier_conflict=True`, and leaves all affected image cells blank. |
| `qxy.import_cells()` | Imports measurements independently while threshold conflicts await review. |
| `qxy.threshold_from_classifiers()` | Raises an error for direct conflicting JSON definitions; otherwise saves the applied values to `thresholds/classifier_thresholds.tsv`. |
| `qxy.threshold_from_table()` | Raises an error for conflict-marked rows with any missing image value. |
| Resolved table | One numeric value per image is applied by exact match to `adata.obs["Image"]`. Values may differ between images. |
| Reapplication | Rerunning thresholding removes old threshold-dependent cell-type outputs; rerun `qxy.celltype()` to regenerate labels. |

    qxy.threshold_from_table(adata, "reviewed_thresholds.tsv")
    qxy.celltype(adata, "celltype_logic.yaml")

## TMA CoreID Commands

**CoreID assignment:** QXYCell recognizes only QuPath's exact measurement-table `TMA Core` column. `qxy.import_cells()` preserves it and automatically creates categorical `CoreID` values.

| Goal | Command | Notes |
|----|----|----|
| Assign CoreID from measurement metadata | `qxy.assign_core_ids_from_measurements(adata)` | Creates `CoreID` from `TMA Core`. This happens automatically inside `qxy.import_cells()` when that source column exists. |
| Check measurement-derived CoreID assignment | `adata.obs["CoreID"].value_counts()` | Missing cells are excluded automatically by pandas counts/grouping. |
| Filter to cells with CoreID | `adata[adata.obs["CoreID"].notna()]` | Excludes cells without measurement-derived core IDs. |

**GeoJSON matching rule:** GeoJSON files are matched to images by comparing the filename stem (stripped of `.ome` and extension) against the `Image` column value stem. Only files whose stem matches a loaded image contribute annotations or cell polygons. GeoJSON never creates CoreID values.

## CN Analysis

CN analysis characterises each cell by the composition of its local spatial neighbourhood. `qxy.cn_knn` builds a per-cell composition profile, `qxy.cn_kmeans` clusters those profiles into discrete CNs, and `qxy.cn_name` assigns deterministic short labels from the composition. All three require **scikit-learn**, which is included in the core dependencies.

| Function | Purpose | Typical command | Key output |
|----|----|----|----|
| `qxy.cn_knn` | For each cell, finds its *k* nearest spatial neighbours within the same image and computes the frequency of each cell type among those neighbours. KNN graphs are built independently per image to prevent neighbours bleeding across tissue sections. | `qxy.cn_knn(adata, k=12)` | Adds `adata.obsm["cn_profile"]` (float32, shape *n_cells × n_cell_types*) and `adata.uns["cn"]`. |
| `qxy.cn_kmeans` | Clusters `adata.obsm["cn_profile"]` into *n_cn* CNs using MiniBatchKMeans. Cluster IDs are assigned as `"N0"` … `"N{n-1}"`. | `qxy.cn_kmeans(adata, n_cn=12)` | Adds `adata.obs["cn"]` as a `pd.Categorical` and updates `adata.uns["cn"]`. |
| `qxy.cn_name` | Replaces cluster IDs (`N0`, `N1`, …) with deterministic short labels derived from the CN composition profile. Rules applied in order: top cell type ≥ 50 % → `<top> hi`; ≥ 35 % → `<top> mix`; top two together ≥ 55 % → `<top1>/<top2>`; otherwise → `<top1>/<top2> mix`. Duplicate labels are disambiguated automatically. An optional `compaction` dict maps long cell type names to shorter display labels before naming. | `qxy.cn_name(adata)` | Updates `adata.obs["cn"]` in place. Saves `cn_labels.csv` and stores `label_map` in `adata.uns["cn"]`. |

### CN Parameters

| Parameter | Function | Default | Notes |
|----|----|----|----|
| `k` | `cn_knn` | `12` | Number of nearest neighbours per cell. Self is excluded. |
| `cell_type_col` | `cn_knn` | `"celltype"` | Column in `adata.obs` containing cell type labels, as assigned by `qxy.celltype()`. |
| `image_col` | `cn_knn` | `"Image"` | Column used to separate cells by tissue section for per-image KNN. |
| `n_cn` | `cn_kmeans` | `12` | Number of CNs to cluster into. `n_clusters` is still accepted as a legacy alias. |
| `key` | `cn_kmeans`, `cn_name` | `"cn"` | Column in `adata.obs` to read/write CN labels. |
| `random_state` | `cn_kmeans` | `0` | Random seed for reproducibility. |
| `compaction` | `cn_name` | `None` | Optional dict mapping original cell type names to shorter display labels, e.g. `{"CD8+PD1+LAG3+": "PD1 LAG3 CD8"}`. Applied before naming. |

### CN Workflow

    # After qxy.celltype(adata) has assigned celltype labels:
    qxy.cn_knn(adata, k=12)
    qxy.cn_kmeans(adata, n_cn=12)
    qxy.cn_name(adata)

    # Optional: compact long cell type names before naming
    compaction = {"CD8+PD1+LAG3+": "PD1 LAG3 CD8", "Mac_MHCII_CD11c-_LY6C-": "MHCII Mac"}
    qxy.cn_name(adata, compaction=compaction)

    # Inspect CN labels
    adata.obs["cn"].value_counts()

    # Spatial plot coloured by CN
    qxy.plot_spatial(adata, sample_col="ImageID", category_col="cn")

    # Cross-tabulate CNs with cell types
    import pandas as pd
    pd.crosstab(adata.obs["cn"], adata.obs["celltype"], normalize="index")

## Plotting Commands

| Goal | Command | Output |
|----|----|----|
| Stacked bar by sample | `qxy.plot_stacked_bar(adata, sample_col="ImageID")` | PNG, PDF, and CSV frequency tables. |
| Stacked bar by group | `qxy.plot_stacked_bar(adata, sample_col="ImageID", group_col="group")` | Group-level mean frequencies. |
| Limit to selected celltypes | `qxy.plot_stacked_bar(adata, celltypes=["T_cell", "Macrophage"])` | Only selected categories plotted. |
| Stacked bar — journal double-column width | `qxy.plot_stacked_bar(adata, width="double")` | 180 mm minimum figure width. |
| Stacked bar — wider bars | `qxy.plot_stacked_bar(adata, bar_width_mm=20)` | 20 mm per bar column instead of the 15 mm default. |
| Spatial plot with automatic grouping | `qxy.plot_spatial(adata)` | Uses usable `Sample` labels when available, otherwise falls back to `Image`. |
| Spatial plot forced by image | `qxy.plot_spatial(adata, sample_col="Image")` | One spatial figure per QuPath image, even when `Sample` labels exist. |
| Combined spatial panel | `qxy.plot_spatial(adata, sample_col="ImageID", combined=True, max_cols=2)` | Individual plots plus combined PNG/PDF. |
| Selected samples | `qxy.plot_spatial(adata, sample_col="ImageID", samples=["sample_A", "sample_B"])` | Restricts plots to listed samples. |
| Exclude missing sample labels | `qxy.plot_spatial(adata, sample_col="Sample")` | Cells where `Sample` is NaN are excluded by default. Pass `include_missing_samples=True` to include a `nan` panel. |
| Full-data underlay for a filtered analysis | `qxy.plot_spatial(adata_cn, underlay_adata=adata, sample_col="Sample", category_col="cn")` | Grey cells and panel bounds come from `adata`; coloured CN cells come from `adata_cn`. |
| Smaller points | `qxy.plot_spatial(adata, point_size=0.3, underlay_size=0.1)` | Better for dense cell maps. |
| Palette control | `qxy.plot_spatial(adata, palette="tab20")` | Uses a matplotlib palette. |
| Cell boundary polygons | `qxy.plot_cell_boundaries(adata, sample_col="ImageID", label_celltypes="Tumor", save_pdf=False)` | Draws WKT cell polygons instead of centroid dots. More memory-intensive than dot plots. |
| Annotation polygon QC | `qxy.plot_annotation_polygons(adata)` | Reloads annotation geometry from the original project, draws a low-resolution cell-density underlay, uses boundary-only polygons by default (`fill=False`), and saves PNG only. |
| Marker positivity heatmap | `qxy.plot_marker_positivity_heatmap(adata)` | Fraction positive per cell type × marker. PDF plus CSV by default; optional PNG/SVG/TIFF. Clustered rows + cols. Batlow colormap (cividis fallback). Colorbar label: `f`. |
| Marker intensity heatmap | `qxy.plot_marker_intensity_heatmap(adata)` | Z-score mean intensity per marker column. Coolwarm colormap, centred 0, ±3. Colorbar label: `z (↓)`. |
| Both marker heatmaps | Call both explicit marker heatmap functions | Produces the two distinct quantities explicitly; the legacy combined entry point is retained only for old code. |
| Marker heatmap — journal single column | `qxy.plot_marker_positivity_heatmap(adata, width="single")` | Fixed 90 mm figure width. |
| Marker heatmap — journal double column | `qxy.plot_marker_intensity_heatmap(adata, width="double")` | Fixed 180 mm figure width. |
| Marker heatmap — selected markers, no cluster | `qxy.plot_marker_intensity_heatmap(adata, markers=["CD45", "CD4", "CD8"], cluster_cols=False)` | Columns in specified order. |
| CN abundance heatmap (by sample) | `qxy.plot_cn_heatmap(adata)` | CN fraction per `Sample`, falling back to `Image` when usable sample labels are unavailable (columns sum to 1). Clustered rows + cols. Batlow colormap. Colorbar label: `f (↓)`. PDF plus CSV by default; optional PNG/SVG/TIFF. |
| CN heatmap with condition strip | `qxy.plot_cn_heatmap(adata, condition_col="group")` | Glasby colour strip above samples grouping by condition, with legend. |
| CN heatmap — both normalisations | `qxy.plot_cn_heatmap(adata, normalize="both")` | Generates sample-normalised and CN-normalised plots. |

### Stacked Bar Parameters

Figure width is determined by `bar_width_mm × n_bars`, with `width="single"` or `"double"` acting as a minimum (same convention as the heatmaps). Colours are shared with `plot_spatial` and the marker heatmap functions via the glasbey palette cached in `adata.uns["qxycell"]["palettes"]`.

| Parameter | Default | Notes |
|----|----|----|
| `width` | `"single"` | `"single"` — 90 mm minimum total figure width. `"double"` — 180 mm minimum. Figure grows beyond the minimum when there are many bars. |
| `bar_width_mm` | `15.0` | Width of each bar column in mm. Total plot width = `bar_width_mm × n_bars` + legend, subject to the `width` minimum. |
| `height_mm` | `72.0` | Total figure height in mm. |
| `legend_width_mm` | `38.0` | Width of the legend panel in mm. |
| `denominator` | `"all_cells"` | Frequency denominator. `"all_cells"` — divide by all cells in the sample (including unplotted categories). `"plotted_cells"` — divide by only the plotted categories. |
| `dpi` | `600` | Resolution for PNG output. PDF is vector. |

### Spatial and Annotation-Polygon Parameters

| Function | Parameter | Default | Meaning |
|----|----|----|----|
| `plot_spatial` | `underlay_adata` | `None` | Use another AnnData for grey cells and panel bounds; the primary `adata` supplies coloured categories. |
| `plot_spatial` | `sample_col` | `None` | Automatically prefer usable `Sample` labels shared by the plotted and underlay objects, otherwise use `Image`. Pass a column name to force it. |
| `plot_spatial` | `include_missing_samples` | `False` | Exclude cells where `adata.obs[sample_col]` is NaN. Set true to include a `nan` panel. |
| `plot_spatial` | `save_png` / `save_pdf` | `True` / `False` | PNG-only by default; enable either format or both. |
| `plot_annotation_polygons` | `project_dir` | `None` | Falls back to `adata.uns["qxycell"]["project_dir"]`. |
| `plot_annotation_polygons` | `images` | `None` | Plot every image unless a list is supplied. |
| `plot_annotation_polygons` | `cell_underlay` | `True` | Draw a coarse 2D cell-density raster beneath annotations. |
| `plot_annotation_polygons` | `underlay_bins` | `384` | Density-raster resolution. |
| `plot_annotation_polygons` | `underlay_cmap` / `underlay_alpha` | `"Greys"` / `0.45` | Cell-density appearance. |
| `plot_annotation_polygons` | `fill` / `fill_alpha` | `False` / `0.20` | Boundary-only by default; enable translucent fills explicitly. |
| `plot_annotation_polygons` | `boundary_linewidth` | `1.0` | Annotation outline width. |
| `plot_annotation_polygons` | `flip_y` | `True` | Match image-viewer orientation. |
| `plot_annotation_polygons` | `figsize` / `dpi` | `(10, 10)` / `300` | PNG size and resolution. |

### Heatmap Parameters

The heatmap functions save **one PDF plus the plotted CSV matrix by default**. PNG, SVG, and TIFF can be enabled individually. Rows and columns are independently reordered by hierarchical clustering (average linkage, Euclidean distance) — no dendrogram is drawn. When `row_strip=True`, a glasbey colour strip is drawn to the left of the heatmap using the same per-category palette as `plot_spatial` and `plot_stacked_bar`, so colours are consistent across all three plot types. Text colour inside annotated cells is chosen for luminance contrast against the background.

Colour consistency: the first call to `plot_spatial`, `plot_stacked_bar`, or either explicit marker heatmap function (with `row_strip=True`) generates a glasbey palette for all categories and caches it in `adata.uns["qxycell"]["palettes"]`. Subsequent calls re-use the same colours. To regenerate (e.g. after adding new cell types): `adata.uns["qxycell"]["palettes"].pop("celltype")`.

| Parameter | Function | Default | Notes |
|----|----|----|----|
| `category_col` | Both marker heatmap functions | `"celltype"` | Column in `adata.obs` used as rows (one row per category). Pass `"cn"` to show markers by CN instead. |
| `markers` | Both marker heatmap functions | `None` | List of marker names to include (column order). By default, only markers actually thresholded into matching `<marker>_pos` columns are used. An explicit list can intentionally include raw, unthresholded intensities. |
| Plotted quantity | `plot_marker_positivity_heatmap` | Fraction positive | Mean of binary `_pos` columns: positive-cell count divided by category cell count, scale 0–1. |
| Plotted quantity | `plot_marker_intensity_heatmap` | Z-scored mean intensity | Mean `adata.X` intensity per category and marker, then Z-scored separately down each marker column. This is not median intensity or a cell count. |
| `cn_col` | `plot_cn_heatmap` | `"cn"` | Column in `adata.obs` containing CN labels (rows). |
| `sample_col` | `plot_cn_heatmap` | `None` | Automatically use usable `Sample` labels, otherwise `Image`. Pass a column name to force it. |
| `condition_col` | `plot_cn_heatmap` | `None` | Optional column mapping samples to experimental conditions. When provided, a glasbey-coloured strip is drawn above the heatmap grouping samples by condition, with a legend. |
| `normalize` | `plot_cn_heatmap` | `"sample"` | `"sample"` — columns sum to 1 (CN fraction per sample); colorbar label `f (↓)` — compare CNs within a sample column. `"cn"` — rows sum to 1 (sample distribution per CN); colorbar label `f (→)` — compare samples within a CN row. `"both"` — generates both plots in one call. Colormap defaults to `batlow` (cividis fallback). |
| `cluster_rows` | both | `True` | Reorder rows by hierarchical clustering. Set `False` to preserve input order. |
| `cluster_cols` | both | `True` | Reorder columns by hierarchical clustering. Set `False` to preserve input order. |
| `width` | both | `"single"` | `"single"` — 90 mm, one journal column (default). `"double"` — 180 mm, two journal columns / full page. `"auto"` — tile-based (0.25 cm × 0.25 cm per cell); produces very small figures for typical marker/CN counts — use only when exact tile size matters. |
| `cmap` | both | `None` | Override the default colormap. Sequential plots (positivity, CN) default to `batlow` if `cmcrameri` is installed, otherwise `cividis`. Intensity defaults to `coolwarm`. Pass any matplotlib colormap name *or* a Crameri short alias (see table below). |
| `annotate` | both | `False` | Write values inside each cell. Text colour (black/white) is chosen automatically for luminance contrast. Default off — pass `annotate=True` to enable. |
| `row_strip` | both | `False` | Draw a glasbey colour strip to the left of the heatmap, one colour per row category. Uses the shared palette cached in `adata.uns["qxycell"]["palettes"]` so colours match `plot_spatial` and `plot_stacked_bar`. |
| `dpi` | both | `600` | Resolution for the TIFF output only. PDF and SVG are vector formats and ignore this parameter. |

### Crameri Colormap Aliases

QXYCell resolves a set of short Crameri colormap names automatically. If `cmcrameri` is installed the full scientific colormap is used; otherwise a comparable matplotlib colormap is substituted silently. Pass these names as `cmap=` to any heatmap function.

| Alias | Type | Crameri name | Fallback | Good for |
|----|----|----|----|----|
| `"batlow"` | Sequential | `cmc.batlow` | `cividis` | Positivity, CN abundance (default) |
| `"batlowS"` | Sequential (cyclic) | `cmc.batlowS` | `tab20` | Categorical / cyclic data |
| `"roma"` | Diverging | `cmc.roma` | `coolwarm` | Z-score intensity (blue→white→red) |
| `"vik"` | Diverging | `cmc.vik` | `coolwarm` | Z-score intensity (cooler tones) |

Install `cmcrameri` to access the full scientific colormaps: `pip install cmcrameri`. Without it, the fallback column is used automatically — no code changes needed.

## Useful AnnData Commands

| Command | Purpose |
|----|----|
| `adata` | Display AnnData summary. |
| `adata.n_obs` | Number of cells. |
| `adata.n_vars` | Number of marker features in `adata.X`. |
| `adata.obs.head()` | Preview cell metadata. |
| `adata.obs.columns.tolist()` | List all cell metadata columns. |
| `adata.var` | View marker metadata. |
| `adata.var_names.tolist()` | List imported markers. |
| `adata.X` | Marker measurement matrix. |
| `adata.obsm["spatial"]` | Spatial coordinates. |
| `adata.uns` | Run metadata and QXYCell summaries. |
| `adata.uns["qxycell"]["output_dir"]` | Current QXYCell output folder. |
| `adata.uns["qxycell"]["h5ad_path"]` | Current H5AD save path. |
| `adata.obs["CoreID"].value_counts()` | Cells per measurement-derived core ID (NaN cells excluded automatically). |
| `adata.obs["CoreID"].cat.categories` | List assigned CoreID names. |
| `adata[adata.obs["CoreID"].notna()]` | Subset to cells with measurement-derived CoreID. |
| `adata.obs["Sample"].value_counts()` | Cells per sample (NaN unassigned cells excluded automatically). |
| `adata.obs["celltype"].value_counts()` | Cell type counts (after `qxy.celltype()`). |
| `adata[adata.obs["celltype"] == "T_cell"].copy()` | Subset to one cell type. |
| `adata.obs["cn"].value_counts()` | Cells per CN (after `qxy.cn_kmeans()`). |
| `adata.obsm["cn_profile"]` | CN composition profile matrix (after `qxy.cn_knn()`). |
| `adata.write_h5ad("file.h5ad")` | Manual H5AD save. |
| `import anndata as ad; adata = ad.read_h5ad("file.h5ad")` | Manual H5AD load. |

## Useful Pandas Commands For `adata.obs`

| Command | Purpose |
|----|----|
| `adata.obs["Image"].value_counts()` | Cells per QuPath image. |
| `adata.obs["ImageID"].value_counts()` | Cells per short sample ID. |
| `adata.obs["celltype"].value_counts()` | Cell type counts. |
| `adata.obs.filter(regex="_pos$").mean()` | Marker positivity fractions. |
| `adata.obs.filter(regex="^annotation__").sum()` | Annotation cell counts. |
| `adata.obs[["Image", "ImageID"]].drop_duplicates()` | Check image-to-sample mapping. |
| `adata.obs.groupby("ImageID").size()` | Cells per sample. |
| `adata.obs.groupby(["ImageID", "celltype"]).size()` | Cell types per sample. |
| `pd.read_csv("metadata.tsv", sep="\t")` | Load metadata table. |
| `adata.obs.to_csv("cells_obs.tsv", sep="\t")` | Export cell metadata. |

## Image Name Mapping

    image_name_map = {
        "long_original_image_name_1.ome.tiff": "sample-1",
        "long_original_image_name_2.ome.tiff": "sample-2",
    }

    image = adata.obs["Image"].astype(str)
    adata.obs["ImageID"] = image.map(image_name_map).fillna(image)

## CLI Commands

| Command | Purpose |
|----|----|
| `qxycell check /path/to/qupath_project` | Validate a QuPath project folder from the terminal. |
| `qxycell import-cells /path/to/qupath_project` | Import the QuPath cell measurement table into AnnData. |
| `qxycell import-cells /path/to/qupath_project --out my_output` | Import cells into an explicit output folder. |

## Output Locations

| Output | Default location |
|----|----|
| Run folder | `<project-name>_run_YYMMDD_HHMM/` beside the project folder |
| Check folder | `<project-name>_check_YYMMDD_HHMM/` beside the project folder |
| H5AD | `<run-folder>/h5ad/qxycell.h5ad` |
| Cell metadata table | `<run-folder>/tables/cells_obs.csv` |
| Marker metadata table | `<run-folder>/tables/markers_var.csv` |
| Annotation assignment audit | `<run-folder>/tables/annotation_assignments.csv` |
| Celltype prompt and summaries | `<run-folder>/celltype/` |
| Dataset summary | `<run-folder>/dataset_summary/dataset_summary.html` |
| Plots | `<run-folder>/plots/` |
