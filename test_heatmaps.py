"""
test_heatmaps.py — QuXYCell heatmap smoke tests
Run from the repo root in your quxycell-dev environment:

    python test_heatmaps.py

Or cell-by-cell in a VSCode Jupyter session (# %% cells below).
"""

# %% ── 0. Imports ────────────────────────────────────────────────────────────

import quxycell as qxy

# Load the most recent run from the TseIvanhoe SLI#3 analysis
QUPATH_DIR = "/Users/json/Desktop/1_Projects/260527_TseIvanhoe/Comet SLI#3 Analysis"
H5AD_PATH  = "/Users/json/Desktop/1_Projects/260527_TseIvanhoe/Comet SLI#3 Analysis/qxy_outputs_260528-0917/run/h5ad/quxycell_260528-0917.h5ad"

adata = qxy.load(H5AD_PATH)
print(adata)
print("\nobs columns:", list(adata.obs.columns))
print("var_names:  ", list(adata.var_names))


# %% ── 1. Marker positivity heatmap (default) ───────────────────────────────
# Rows = cell types, Cols = markers
# Colormap: cividis (0 → 1 fraction positive)
# Rows + cols clustered. Left colour strip per cell type (glasbey palette).
# Saves: PDF, SVG, TIFF @ 300 dpi, CSV

qxy.plot_marker_heatmap(adata)


# %% ── 2. Marker intensity heatmap ──────────────────────────────────────────
# Colormap: coolwarm, centred at 0, clipped ±3

qxy.plot_marker_heatmap(adata, values="intensity")


# %% ── 3. Both modes in one call ─────────────────────────────────────────────

qxy.plot_marker_heatmap(adata, values="both")


# %% ── 4. Journal single-column width (90 mm) ───────────────────────────────
# This is now the default — both lines are equivalent:

qxy.plot_marker_heatmap(adata)
qxy.plot_marker_heatmap(adata, width="single")


# %% ── 5. Journal double-column width (180 mm) ──────────────────────────────

qxy.plot_marker_heatmap(adata, width="double")


# %% ── 6. Selected markers only, columns in fixed order (no col clustering) ─

qxy.plot_marker_heatmap(
    adata,
    markers=["CD_45", "CD4", "CD14", "CD56", "PD_1"],
    cluster_cols=False,
)


# %% ── 7. Override colormap ──────────────────────────────────────────────────

qxy.plot_marker_heatmap(adata, cmap="YlOrRd")


# %% ── 8. Turn off annotations (cleaner for many cell types) ─────────────────

qxy.plot_marker_heatmap(adata, annotate=False)


# %% ── 9. Group by CN instead of celltype ────────────────────────────────────
# Useful to see which markers define each CN

qxy.plot_marker_heatmap(adata, category_col="cn")


# %% ── 10. CN abundance heatmap — fraction per sample (default) ──────────────
# Rows = CNs, Cols = samples
# Colormap: cividis. Left CN colour strip. Clustered rows + cols.

qxy.plot_cn_heatmap(adata)


# %% ── 11. CN heatmap — flip normalisation (fraction of each CN per sample) ──

qxy.plot_cn_heatmap(adata, normalize="cn")


# %% ── 12. Both normalisation modes in one call ───────────────────────────────

qxy.plot_cn_heatmap(adata, normalize="both")


# %% ── 13. CN heatmap with a condition strip ─────────────────────────────────
# If your adata.obs has a "group" or "condition" column, a coloured strip
# appears above the heatmap grouping samples by condition, with a legend.
# Replace "group" with whatever column you have.

# qxy.plot_cn_heatmap(adata, condition_col="group")


# %% ── 14. CN heatmap — use short sample names, single-column width ──────────

qxy.plot_cn_heatmap(adata, sample_col="Image", width="single")


# %% ── 15. CN heatmap — double-column width, no annotations ──────────────────

qxy.plot_cn_heatmap(adata, width="double", annotate=False)


# %% ── 16. Override colormap on CN heatmap ───────────────────────────────────

qxy.plot_cn_heatmap(adata, cmap="Purples")


# %% ── 17. Disable clustering (preserve original row/col order) ──────────────

qxy.plot_marker_heatmap(adata, cluster_rows=False, cluster_cols=False)
qxy.plot_cn_heatmap(adata, cluster_rows=False, cluster_cols=False)


# %% ── 18. Save to a specific output directory ───────────────────────────────
# By default plots go into the folder stored in adata.uns["quxycell"]["output_dir"].
# Override with output_dir= to redirect elsewhere.

import pathlib
custom_dir = pathlib.Path(QUPATH_DIR) / "qxy_outputs_260528-0917"

qxy.plot_marker_heatmap(adata, output_dir=custom_dir, show=False)
qxy.plot_cn_heatmap(adata, output_dir=custom_dir, show=False)


# %% ── 19. Suppress inline display (for batch / headless runs) ───────────────

qxy.plot_marker_heatmap(adata, values="both", show=False, verbose=True)
qxy.plot_cn_heatmap(adata, normalize="both", show=False, verbose=True)


# %% ── 20. Check what was saved ───────────────────────────────────────────────

result = qxy.plot_marker_heatmap(adata, show=False, verbose=False)
for mode, paths in result.items():
    print(f"\n{mode}:")
    for p in paths:
        print(f"  {p}")
