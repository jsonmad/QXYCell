"""Run the complete QXYCell workflow as explicit, editable steps."""

from pathlib import Path

import qxycell as qxy


# 1. CONFIGURE PATHS AND ANALYSIS OPTIONS
PROJECT_DIR = Path("/replace/with/your/qupath/project")
OUTPUT_DIR: Path | None = None
THRESHOLD_FILE: Path | None = None
METADATA_FILE: Path | None = None
CELLTYPE_LOGIC: Path | None = None
SAMPLE_COL = "Sample"
METADATA_SAMPLE_COL: str | None = None
RUN_METADATA = False
RUN_CELLTYPING = False
RUN_PLOTS = True
RUN_CN_ANALYSIS = True
REMOVE_IGNORE_CELLS = True
K_NEIGHBORS = 12
N_CELLULAR_NEIGHBORHOODS = 12


def main() -> None:
    """Run each QXYCell stage in order after editing the configuration above."""

    if not PROJECT_DIR.exists():
        raise FileNotFoundError(
            "Edit PROJECT_DIR at the top of this script before running it: "
            f"{PROJECT_DIR}"
        )

    # 2. CHECK THE QUPath EXPORT
    # check() validates the project and writes a report; it does not create AnnData.
    # Options: project_dir (required QuPath export directory); output_dir=None
    # (report folder); count_rows=False (count CSV rows); threshold_file=None
    # (explicit threshold table, otherwise discovery is used).
    report = qxy.check(
        PROJECT_DIR,
        output_dir=OUTPUT_DIR,
        count_rows=False,
        threshold_file=THRESHOLD_FILE,
    )
    if not report.ok:
        raise RuntimeError("QXYCell project validation failed; inspect the check report above.")

    # 3. IMPORT MEASUREMENTS, SPATIAL DATA, AND ANNOTATIONS
    # run() creates AnnData with measurements, spatial coordinates, annotations,
    # polygons, and run metadata. Thresholding and cell typing stay explicit below.
    # Options: project_dir (required export directory); output_dir=None (run
    # folder); fail_on_check_error=True (stop on validation errors);
    # pixel_size_um=0.28 (coordinate scale); threshold_file=None (threshold
    # table); apply_thresholds=False (do not threshold during import);
    # celltype_logic=None (do not type cells during import); verbose=True
    # (print progress).
    adata = qxy.run(
        PROJECT_DIR,
        output_dir=OUTPUT_DIR,
        fail_on_check_error=True,
        pixel_size_um=0.28,
        threshold_file=THRESHOLD_FILE,
        apply_thresholds=False,
        celltype_logic=None,
        verbose=True,
    )

    # 4. APPLY MARKER THRESHOLDS
    # threshold() adds <marker>_pos columns to adata.obs from the threshold table.
    # Options: adata (imported AnnData); project_dir=None (project path, inferred
    # from adata when omitted); threshold_file=None (explicit table);
    # output_dir=None (output folder); image_col="Image" (per-image threshold
    # column); verbose=True (print progress).
    qxy.threshold(
        adata,
        project_dir=PROJECT_DIR,
        threshold_file=THRESHOLD_FILE,
        output_dir=OUTPUT_DIR,
        image_col="Image",
        verbose=True,
    )

    # 5. OPTIONALLY ADD SAMPLE METADATA
    # add_metadata() joins selected metadata columns onto adata.obs.
    # Options: adata (imported AnnData); metadata (required metadata file or
    # table); sample_col="Image" (AnnData sample column);
    # metadata_sample_col=None (matching metadata column); columns=None (all
    # metadata columns); prefix="" (new-column prefix); overwrite=False
    # (preserve existing columns); output_dir=None (summary folder); verbose=True
    # (print progress).
    if RUN_METADATA and METADATA_FILE is not None:
        qxy.add_metadata(
            adata,
            METADATA_FILE,
            sample_col=SAMPLE_COL,
            metadata_sample_col=METADATA_SAMPLE_COL,
            columns=None,
            prefix="",
            overwrite=False,
            output_dir=OUTPUT_DIR,
            verbose=True,
        )
    else:
        print("Skipping metadata: set RUN_METADATA=True and METADATA_FILE to enable it.")

    # 6. REMOVE IGNORE-ANNOTATED CELLS
    # remove_ignore() removes cells flagged by Ignore annotations. With copy=False
    # it mutates adata in place, so its return value is intentionally not assigned.
    # Options: adata (imported AnnData); annotation_prefix="annotation__"
    # (annotation-column prefix); ignore_text="Ignore" (annotation text to
    # remove); copy=False (mutate in place); verbose=True (print progress).
    if REMOVE_IGNORE_CELLS:
        qxy.remove_ignore(
            adata,
            annotation_prefix="annotation__",
            ignore_text="Ignore",
            copy=False,
            verbose=True,
        )
    else:
        print("Skipping Ignore-cell removal: REMOVE_IGNORE_CELLS=False.")

    # 7. OPTIONALLY ASSIGN CELL TYPES
    # celltype() applies ordered rules and writes the celltype column and any
    # feature flags to adata.obs.
    # Options: adata (imported AnnData); logic=None (YAML path or rule mapping);
    # celltype_column="celltype" (destination column); unknown_label="Unknown"
    # (label for unmatched cells); celltype_dir=None (default rule/output folder);
    # verbose=True (print progress).
    if RUN_CELLTYPING and CELLTYPE_LOGIC is not None:
        qxy.celltype(
            adata,
            logic=CELLTYPE_LOGIC,
            celltype_column="celltype",
            unknown_label="Unknown",
            celltype_dir=None,
            verbose=True,
        )
    else:
        print("Skipping cell typing: set RUN_CELLTYPING=True and CELLTYPE_LOGIC to enable it.")

    # 8. WRITE QUALITY-CONTROL SUMMARIES
    # qc() writes tabular and HTML summaries for cells, annotations, markers, and
    # cell types that are present in adata.
    # Options: adata (processed AnnData); sample_col="Image" (sample column);
    # celltype_col="celltype" (cell-type column); annotation_prefix="annotation__"
    # (annotation-column prefix); marker_suffix="_pos" (positivity-column suffix);
    # output_dir=None (QC folder); verbose=True (print progress).
    qxy.qc(
        adata,
        sample_col=SAMPLE_COL,
        celltype_col="celltype",
        annotation_prefix="annotation__",
        marker_suffix="_pos",
        output_dir=OUTPUT_DIR,
        verbose=True,
    )

    # 9. CREATE STANDARD PLOTS
    # These plots use celltype categories, so they require the celltype column.
    if RUN_PLOTS and "celltype" in adata.obs:
        # plot_stacked_bar() plots category frequencies by sample.
        # Options: adata (processed AnnData); category_col="celltype" (stack
        # labels); sample_col="Image" (bar groups); group_col=None (optional
        # aggregation); subset_col=None and subset_value=None (optional filter);
        # samples=None and celltypes=None (optional selections); output_dir=None
        # (plot folder); filename_prefix=None and save_prefix=None (file names);
        # colors=None and palette=None (category colours); denominator="all_cells"
        # (frequency denominator); width="single" (minimum figure width);
        # bar_width_mm=15.0, height_mm=72.0, legend_width_mm=38.0 (dimensions);
        # show_axis_labels=True, x_axis_label=None, y_axis_label=None (axis
        # labels); dpi=600 (PNG resolution); save_png=False and save_pdf=True
        # (PDF output); show=True (display figure); verbose=True (print progress).
        qxy.plot_stacked_bar(
            adata,
            category_col="celltype",
            sample_col=SAMPLE_COL,
            output_dir=OUTPUT_DIR,
            save_png=False,
            save_pdf=True,
            show=False,
            verbose=True,
        )

        # plot_spatial() plots category-coloured cell centroids in each sample.
        # Options: adata (processed AnnData); underlay_adata=None (all-cell
        # underlay source); category_col="celltype" (colour categories);
        # sample_col="Image" (sample column); subset_col=None and
        # subset_value=None (optional filter); samples=None, celltypes=None, and
        # images=None (optional selections); include_missing_samples=False
        # (exclude missing labels); spatial_key=None (coordinate key);
        # output_dir=None (plot folder); filename_prefix=None and save_prefix=None
        # (file names); colors=None and palette=None (category colours);
        # fixed_window_um=None and center_method="bbox" (view framing);
        # point_size=4.0, underlay_size=2.0, underlay_color="#bdbdbd", and
        # underlay_alpha=0.08 (point and underlay style); scale_bar=True,
        # scale_bar_um=1000.0, scale_bar_label="1 mm" (scale bar); flip_y=True
        # (image orientation); figsize=(10.0, 10.0), auto_figsize=False,
        # dpi=600, legend_width=2.2 (layout); combined=False,
        # save_individual=True, save_png=True, save_pdf=False, max_cols=3
        # (PNG output/layout); show=True (display); verbose=True (print progress).
        qxy.plot_spatial(
            adata,
            category_col="celltype",
            sample_col=SAMPLE_COL,
            output_dir=OUTPUT_DIR,
            save_png=True,
            save_pdf=False,
            show=False,
            verbose=True,
        )

        # plot_annotation_polygons() overlays QuPath annotation polygons for QC.
        # Options: adata (processed AnnData); project_dir=None (GeoJSON project
        # directory); image_col="Image" (image column); images=None (selected
        # images); output_dir=None (plot folder); pixel_size_um=None (coordinate
        # scale); colors=None (annotation colours); cell_underlay=True,
        # underlay_bins=384, underlay_cmap="Greys", underlay_alpha=0.45
        # (underlay style); fill=False, fill_alpha=0.2, boundary_linewidth=1.0
        # (polygon style); flip_y=True (image orientation); figsize=(10.0, 10.0),
        # dpi=300 (PNG layout); show=True (display); verbose=True (print progress).
        # This function writes PNG-only annotation-polygon QC output.
        qxy.plot_annotation_polygons(
            adata,
            project_dir=PROJECT_DIR,
            image_col=SAMPLE_COL,
            output_dir=OUTPUT_DIR,
            pixel_size_um=0.28,
            show=False,
            verbose=True,
        )

        # plot_marker_positivity_heatmap() shows per-category marker positivity.
        # Options: adata (processed AnnData); category_col="celltype" (row
        # categories); markers=None (all thresholded markers); cluster_rows=True
        # and cluster_cols=True (hierarchical ordering); width="single" (figure
        # width); cmap=None (default colour map); annotate=False (cell labels);
        # row_strip=False (row-colour strip); dpi=600 (raster resolution);
        # save_png=False, save_pdf=True, save_svg=False, save_tiff=False
        # (PDF output); output_dir=None (plot folder); show=True (display);
        # verbose=True (print progress).
        qxy.plot_marker_positivity_heatmap(
            adata,
            category_col="celltype",
            output_dir=OUTPUT_DIR,
            save_png=False,
            save_pdf=True,
            save_svg=False,
            save_tiff=False,
            show=False,
            verbose=True,
        )

        # plot_marker_intensity_heatmap() shows Z-scored category mean intensity.
        # Options: adata (processed AnnData); category_col="celltype" (row
        # categories); markers=None (all thresholded markers); cluster_rows=True
        # and cluster_cols=True (hierarchical ordering); width="single" (figure
        # width); cmap=None (default colour map); annotate=False (cell labels);
        # row_strip=False (row-colour strip); dpi=600 (raster resolution);
        # save_png=False, save_pdf=True, save_svg=False, save_tiff=False
        # (PDF output); output_dir=None (plot folder); show=True (display);
        # verbose=True (print progress).
        qxy.plot_marker_intensity_heatmap(
            adata,
            category_col="celltype",
            output_dir=OUTPUT_DIR,
            save_png=False,
            save_pdf=True,
            save_svg=False,
            save_tiff=False,
            show=False,
            verbose=True,
        )

        # plot_cell_boundaries() colours cells by their WKT boundary polygons.
        # Options: adata (processed AnnData); category_col="celltype" (colour
        # categories); sample_col="Image" (sample column);
        # polygon_col="cell_polygon_wkt" (WKT column); subset_col=None and
        # subset_value=None (optional filter); samples=None, celltypes=None, and
        # images=None (optional selections); spatial_key=None (coordinate key);
        # output_dir=None (plot folder); filename_prefix=None and save_prefix=None
        # (file names); colors=None and palette=None (category colours);
        # fixed_window_um=None and center_method="bbox" (view framing); fill=True,
        # fill_alpha=0.85, boundary_linewidth=0.08, boundary_color=None
        # (polygon style); underlay=True, underlay_facecolor="#d9d9d9",
        # underlay_edgecolor="#bdbdbd", underlay_alpha=0.25,
        # underlay_linewidth=0.04 (underlay style); scale_bar=True,
        # scale_bar_um=1000.0, scale_bar_label="1 mm" (scale bar); flip_y=True
        # (image orientation); label_celltypes=None, label_max_per_celltype=1,
        # label_offset_um=(150.0, 150.0), label_fontsize=7.0,
        # label_linewidth=0.6, label_color=None (labels); figsize=(10.0, 10.0),
        # auto_figsize=False, dpi=600, legend_width=2.2 (layout);
        # save_individual=True, save_png=True, save_pdf=False (PNG output);
        # show=True (display); verbose=True (print progress).
        if "cell_polygon_wkt" in adata.obs:
            qxy.plot_cell_boundaries(
                adata,
                category_col="celltype",
                sample_col=SAMPLE_COL,
                polygon_col="cell_polygon_wkt",
                output_dir=OUTPUT_DIR,
                save_png=True,
                save_pdf=False,
                show=False,
                verbose=True,
            )
        else:
            print("Skipping cell-boundary plot: cell_polygon_wkt is not available.")
    elif RUN_PLOTS:
        print("Skipping cell-type plots: celltype is not available.")
    else:
        print("Skipping standard plots: RUN_PLOTS=False.")

    # 10. OPTIONALLY ANALYSE CELLULAR NEIGHBOURHOODS
    if RUN_CN_ANALYSIS and "celltype" in adata.obs:
        # cn_knn() stores each cell's same-image neighbour composition in
        # adata.obsm["cn_profile"]. Options: adata (processed AnnData); k=12
        # (neighbours per cell); cell_type_col="celltype" (type column);
        # image_col="Image" (within-image grouping).
        qxy.cn_knn(
            adata,
            k=K_NEIGHBORS,
            cell_type_col="celltype",
            image_col=SAMPLE_COL,
        )

        # cn_kmeans() clusters the profile and adds adata.obs["cn"].
        # Options: adata (profiled AnnData); n_cn=12 (number of CNs);
        # key="cn" (output column); random_state=0 (repeatable clustering);
        # n_clusters=None (legacy alias for n_cn).
        qxy.cn_kmeans(
            adata,
            n_cn=N_CELLULAR_NEIGHBORHOODS,
            key="cn",
            random_state=0,
            n_clusters=None,
        )

        # cn_name() assigns readable composition names to adata.obs["cn"].
        # Options: adata (clustered AnnData); key="cn" (CN column);
        # compaction=None (optional name mapping); output_dir=None (summary
        # folder); verbose=True (print progress).
        qxy.cn_name(
            adata,
            key="cn",
            compaction=None,
            output_dir=OUTPUT_DIR,
            verbose=True,
        )

        # plot_spatial() can also render CN labels. Its full options are the
        # same as the Stage 9 spatial call; this call uses category_col="cn" and
        # writes individual PNGs without opening figures.
        qxy.plot_spatial(
            adata,
            category_col="cn",
            sample_col=SAMPLE_COL,
            output_dir=OUTPUT_DIR,
            save_png=True,
            save_pdf=False,
            show=False,
            verbose=True,
        )

        # plot_cn_heatmap() shows CN composition across samples.
        # Options: adata (CN-labelled AnnData); cn_col="cn" (CN column);
        # category_col=None (alias for cn_col); sample_col="Image" (sample
        # column); include_missing_samples=False (exclude missing labels);
        # condition_col=None (optional sample condition); normalize="sample"
        # (normalisation); cluster_rows=True and cluster_cols=True (ordering);
        # width="single" (figure width); cmap=None (default colour map);
        # annotate=False (cell labels); row_strip=False (row-colour strip);
        # dpi=600 (raster resolution); save_png=False, save_pdf=True,
        # save_svg=False, save_tiff=False (PDF output); output_dir=None (plot
        # folder); show=True (display); verbose=True (print progress).
        qxy.plot_cn_heatmap(
            adata,
            cn_col="cn",
            category_col=None,
            sample_col=SAMPLE_COL,
            include_missing_samples=False,
            condition_col=None,
            normalize="sample",
            cluster_rows=True,
            cluster_cols=True,
            width="single",
            cmap=None,
            annotate=False,
            row_strip=False,
            dpi=600,
            save_png=False,
            save_pdf=True,
            save_svg=False,
            save_tiff=False,
            output_dir=OUTPUT_DIR,
            show=False,
            verbose=True,
        )
    elif RUN_CN_ANALYSIS:
        print("Skipping cellular-neighbourhood analysis: celltype is not available.")
    else:
        print("Skipping cellular-neighbourhood analysis: RUN_CN_ANALYSIS=False.")

    # 11. SAVE THE FINAL AnnData OBJECT
    # save() writes the final H5AD after all enabled optional stages.
    # Options: adata (processed AnnData); path=None (explicit H5AD path);
    # output_dir=None (default output folder); verbose=True (print progress).
    h5ad_path = qxy.save(adata, path=None, output_dir=OUTPUT_DIR, verbose=True)
    print(f"Final QXYCell H5AD: {h5ad_path}")


if __name__ == "__main__":
    main()
