"""Optional step 6: plot assigned cell types in spatial coordinates."""

import qxycell as qxy

from config import OUTPUT_DIR


def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.plot_spatial(
        adata,
        underlay_adata=None,
        category_col="celltype",
        sample_col="Image",
        subset_col=None,
        subset_value=None,
        samples=None,
        celltypes=None,
        images=None,
        include_missing_samples=False,
        spatial_key=None,
        output_dir=None,
        filename_prefix=None,
        save_prefix=None,
        colors=None,
        palette=None,
        fixed_window_um=None,
        center_method="bbox",
        point_size=4.0,
        underlay_size=2.0,
        underlay_color="#bdbdbd",
        underlay_alpha=0.08,
        scale_bar=True,
        scale_bar_um=1000.0,
        scale_bar_label="1 mm",
        flip_y=True,
        figsize=(10.0, 10.0),
        auto_figsize=False,
        dpi=600,
        legend_width=2.2,
        combined=False,
        save_individual=True,
        save_png=True,
        save_pdf=False,
        max_cols=3,
        show=False,  # Library default is True; False suits command-line runs.
        verbose=True,
    )


if __name__ == "__main__":
    main()
