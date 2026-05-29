"""Generic plotting helpers for QXYCell AnnData objects."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Iterable

from qxycell.paths import resolve_output_dir


def _resolve_plot_dir(adata, output_dir: str | Path | None) -> Path:
    if output_dir is not None:
        return Path(output_dir).expanduser().resolve()
    return resolve_output_dir(adata=adata) / "plots"


def _require_plotting():
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mtick
        import numpy as np
        import pandas as pd
        from matplotlib.lines import Line2D
        from matplotlib.colors import hsv_to_rgb, to_hex
    except ImportError as exc:
        raise ImportError(
            "QXYCell plotting requires plot dependencies. Install with "
            "`pip install -e '.[plot]'` or `pip install qxycell[plot]`."
        ) from exc
    return plt, mtick, np, pd, Line2D, hsv_to_rgb, to_hex


def _safe_name(value: object) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(value)).strip("_")
    return safe or "plot"


def _fallback_color(label: object, hsv_to_rgb, to_hex) -> str:
    digest = hashlib.sha1(str(label).encode("utf-8")).hexdigest()
    hue = (int(digest[:8], 16) % 360) / 360.0
    return to_hex(hsv_to_rgb((hue, 0.72, 0.90)))


# Default matplotlib colormap to use per category column name.
# Columns not listed here fall back to the glasbey palette.
_COLUMN_CMAP_DEFAULTS: "dict[str, str]" = {
    "cn": "tab20",
}


def _get_category_palette(adata, category_col: str) -> "dict[str, str]":
    """Return a stable ``{label: hex_colour}`` mapping for every category in
    *category_col*, cached in ``adata.uns["qxycell"]["palettes"]``.

    The palette generator is chosen by ``_COLUMN_CMAP_DEFAULTS``:
    - ``"cn"`` columns use ``tab20``.
    - All other columns use the glasbey palette.

    Caching ensures the same colour is used across ``plot_spatial``,
    ``plot_stacked_bar``, and the heatmap row-strip.
    Call ``adata.uns["qxycell"]["palettes"].pop(category_col)`` to regenerate.
    """
    import matplotlib.pyplot as _plt
    import matplotlib.colors as _mcolors

    uns = adata.uns.setdefault("qxycell", {})
    palettes = uns.setdefault("palettes", {})
    if category_col in palettes:
        return dict(palettes[category_col])

    all_cats = sorted(adata.obs[category_col].astype(str).unique().tolist())
    cmap_name = _COLUMN_CMAP_DEFAULTS.get(category_col)

    palette: dict[str, str] = {}
    if cmap_name is not None:
        cmap = _plt.get_cmap(cmap_name)
        n = len(all_cats)
        for i, cat in enumerate(all_cats):
            palette[cat] = _mcolors.to_hex(cmap(i % cmap.N if hasattr(cmap, "N") else i / max(n - 1, 1)))
    else:
        raw_colors = _glasbey_colors(len(all_cats))
        for cat, color in zip(all_cats, raw_colors):
            palette[cat] = color if isinstance(color, str) else _mcolors.to_hex(color)

    palettes[category_col] = palette
    return dict(palette)


def _color_map(
    labels: Iterable[object],
    hsv_to_rgb,
    to_hex,
    colors: dict[str, str] | None = None,
    palette: dict[str, str] | list[str] | tuple[str, ...] | str | None = None,
):
    colors = colors or {}
    if isinstance(palette, dict):
        colors = {**palette, **colors}
    palette_values = None
    if isinstance(palette, (list, tuple)):
        palette_values = list(palette)
    elif isinstance(palette, str):
        try:
            import matplotlib.pyplot as plt

            cmap = plt.get_cmap(palette)
            labels_list = list(labels)
            denominator = max(len(labels_list) - 1, 1)
            palette_values = [to_hex(cmap(index / denominator)) for index in range(len(labels_list))]
            labels = labels_list
        except Exception:
            palette_values = None

    out = {}
    for index, label in enumerate(labels):
        key = str(label)
        if key in colors:
            out[key] = colors[key]
        elif palette_values:
            out[key] = palette_values[index % len(palette_values)]
        else:
            out[key] = _fallback_color(key, hsv_to_rgb, to_hex)
    return out


def _stack_ymax(df, np, pad: float = 1.08) -> float:
    if df.empty:
        return 0.05
    max_stack = float(df.sum(axis=1).max())
    if not np.isfinite(max_stack) or max_stack <= 0:
        return 0.05
    target = max_stack * pad
    if target <= 0.15:
        step = 0.01
    elif target <= 0.40:
        step = 0.02
    else:
        step = 0.05
    return float(np.ceil(target / step) * step)


def _resolve_spatial_key(adata, spatial_key: str | None) -> str:
    if spatial_key is not None:
        if spatial_key not in adata.obsm:
            raise KeyError(f"Spatial key not found in adata.obsm: {spatial_key}")
        return spatial_key
    if "spatial" in adata.obsm:
        return "spatial"
    if "X_spatial" in adata.obsm:
        return "X_spatial"
    raise KeyError("No spatial coordinates found in adata.obsm. Expected 'spatial' or 'X_spatial'.")


def _add_scale_bar(ax, x_lim, y_lim, *, length_um: float, label: str):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    width = x_max - x_min
    height = y_max - y_min
    if width <= 0 or height <= 0 or length_um <= 0:
        return

    x0 = x_min + width * 0.05
    y0 = y_min + height * 0.05
    x1 = x0 + float(length_um)
    ax.plot([x0, x1], [y0, y0], color="black", linewidth=2.0, solid_capstyle="butt", zorder=5)
    ax.text(
        x0,
        y0 + height * 0.018,
        label,
        color="black",
        fontsize=8,
        ha="left",
        va="bottom",
        zorder=5,
    )


def _prepare_obs(
    adata,
    *,
    category_col: str,
    sample_col: str,
    group_col: str | None = None,
    subset_col: str | None = None,
    subset_value: str | None = None,
    samples: Iterable[str] | None = None,
    categories: Iterable[str] | None = None,
    exclude_categories: Iterable[str] = ("Unknown", "Negative"),
):
    obs = adata.obs.copy()
    required = [category_col, sample_col]
    if group_col is not None:
        required.append(group_col)
    if subset_col is not None:
        required.append(subset_col)
    missing = [column for column in required if column not in obs.columns]
    if missing:
        raise KeyError(f"Missing required adata.obs column(s): {missing}")

    obs = obs[obs[category_col].notna()].copy()
    obs[category_col] = obs[category_col].astype(str).str.strip()
    obs[sample_col] = obs[sample_col].astype(str).str.strip()

    exclude = {str(value).strip().lower() for value in exclude_categories}
    if exclude:
        obs = obs[~obs[category_col].str.lower().isin(exclude)].copy()

    if samples is not None:
        sample_values = {str(sample) for sample in samples}
        obs = obs[obs[sample_col].astype(str).isin(sample_values)].copy()

    if categories is not None:
        category_values = {str(category) for category in categories}
        obs = obs[obs[category_col].astype(str).isin(category_values)].copy()

    if subset_col is not None and subset_value is not None:
        obs = obs[obs[subset_col].astype(str).str.strip() == str(subset_value)].copy()

    return obs


def plot_stacked_bar(
    adata,
    *,
    category_col: str = "celltype",
    sample_col: str = "Image",
    group_col: str | None = None,
    subset_col: str | None = None,
    subset_value: str | None = None,
    samples: Iterable[str] | None = None,
    celltypes: Iterable[str] | None = None,
    output_dir: str | Path | None = None,
    filename_prefix: str | None = None,
    save_prefix: str | None = None,
    colors: dict[str, str] | None = None,
    palette: dict[str, str] | list[str] | tuple[str, ...] | str | None = None,
    denominator: str = "all_cells",
    width: str = "single",
    bar_width_mm: float = 15.0,
    height_mm: float = 72.0,
    legend_width_mm: float = 38.0,
    dpi: int = 600,
    show: bool = True,
    verbose: bool = True,
) -> dict[str, Path]:
    """Create a generic stacked bar plot from AnnData cell annotations.

    Frequencies are computed per sample column, with optional aggregation by
    ``group_col``. Set ``sample_col`` to use a shortened image label column such
    as ``ImageID`` instead of the default QuPath ``Image`` column. By default,
    each category count is divided by all cells in the sample.

    Figure width is determined by ``bar_width_mm`` × number of bars, with
    ``width="single"`` (90 mm) or ``"double"`` (180 mm) acting as a minimum.
    """

    plt, mtick, np, pd, Line2D, hsv_to_rgb, to_hex = _require_plotting()
    output_dir = _resolve_plot_dir(adata, output_dir) / "stacked_bar"
    output_dir.mkdir(parents=True, exist_ok=True)

    obs_plot = _prepare_obs(
        adata,
        category_col=category_col,
        sample_col=sample_col,
        group_col=group_col,
        subset_col=subset_col,
        subset_value=subset_value,
        samples=samples,
        categories=celltypes,
    )
    if obs_plot.empty:
        raise ValueError("No cells available for stacked bar plotting after filtering.")

    obs_all = adata.obs.copy()
    obs_all[sample_col] = obs_all[sample_col].astype(str).str.strip()
    if samples is not None:
        sample_values = {str(sample) for sample in samples}
        obs_all = obs_all[obs_all[sample_col].isin(sample_values)].copy()
    if denominator == "all_cells":
        sample_totals = obs_all.groupby(sample_col, observed=True).size()
    elif denominator == "plotted_cells":
        sample_totals = obs_plot.groupby(sample_col, observed=True).size()
    else:
        raise ValueError("denominator must be 'all_cells' or 'plotted_cells'.")

    counts = (
        obs_plot.groupby([sample_col, category_col], observed=True)
        .size()
        .rename("count")
        .reset_index()
    )
    counts["frequency"] = counts["count"] / counts[sample_col].map(sample_totals)
    freq_sample = (
        counts.pivot(index=sample_col, columns=category_col, values="frequency")
        .fillna(0.0)
    )
    if celltypes is not None:
        category_order = [str(celltype) for celltype in celltypes if str(celltype) in freq_sample.columns]
    else:
        category_order = freq_sample.mean(axis=0).sort_values(ascending=False).index.tolist()
    freq_sample = freq_sample[category_order]

    if group_col is not None:
        sample_groups = (
            obs_plot[[sample_col, group_col]]
            .drop_duplicates()
            .set_index(sample_col)
        )
        freq_plot = (
            freq_sample.join(sample_groups, how="left")
            .groupby(group_col, observed=True)
            .mean(numeric_only=True)
        )
        freq_plot = freq_plot[category_order]
        x_label = group_col
    else:
        freq_plot = freq_sample
        x_label = sample_col

    safe_parts = ["stacked_bar"]
    if subset_col and subset_value:
        safe_parts.append(_safe_name(subset_value))
    if group_col:
        safe_parts.append(f"by_{_safe_name(group_col)}")
    prefix = save_prefix or filename_prefix or "_".join(safe_parts)

    _base_palette = _get_category_palette(adata, category_col)
    color_lookup = _color_map(
        category_order, hsv_to_rgb, to_hex, colors, palette if palette is not None else _base_palette
    )
    stack_colors = [color_lookup[str(label)] for label in category_order]

    # ── Figure sizing (mm → inches) ──────────────────────────────────────────
    _mm = 1.0 / 25.4
    n_bars      = len(freq_plot)
    plot_w_in   = n_bars * bar_width_mm * _mm
    legend_w_in = legend_width_mm * _mm
    total_w_in  = plot_w_in + legend_w_in
    if width == "single":
        total_w_in = max(total_w_in, 90.0 * _mm)
    elif width == "double":
        total_w_in = max(total_w_in, 180.0 * _mm)
    # recompute plot width after any minimum bump
    plot_w_in  = total_w_in - legend_w_in
    height_in  = height_mm * _mm

    fig = plt.figure(figsize=(total_w_in, height_in), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[plot_w_in, legend_w_in], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    legend_ax = fig.add_subplot(gs[0, 1])
    legend_ax.axis("off")
    freq_plot.plot(
        kind="bar",
        stacked=True,
        color=stack_colors,
        ax=ax,
        width=0.8,
        edgecolor="none",
        linewidth=0,
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel("Frequency")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, decimals=0))
    ax.set_ylim(0, _stack_ymax(freq_plot, np))
    ax.tick_params(axis="x", rotation=0)
    ax.set_title(subset_value if subset_value else f"{category_col} frequency", pad=35)
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.get_legend()
    if leg is not None:
        leg.remove()
    legend_ax.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(0, 0.5),
        frameon=False,
        fontsize=8,
        borderaxespad=0.0,
    )
    fig_path = output_dir / f"{prefix}.png"
    pdf_path = output_dir / f"{prefix}.pdf"
    fig.savefig(fig_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)

    table_path = output_dir / f"{prefix}_frequency.csv"
    freq_plot.to_csv(table_path)
    per_image_path = output_dir / f"{prefix}_per_image_frequency.csv"
    freq_sample.to_csv(per_image_path)

    if verbose:
        print(f"Saved stacked bar plot: {fig_path}")
        print(f"Saved stacked bar table: {table_path}")

    return {
        "png": fig_path,
        "pdf": pdf_path,
        "table": table_path,
        "per_image_table": per_image_path,
    }


def plot_spatial(
    adata,
    *,
    category_col: str = "celltype",
    sample_col: str = "Image",
    subset_col: str | None = None,
    subset_value: str | None = None,
    samples: Iterable[str] | None = None,
    celltypes: Iterable[str] | None = None,
    images: Iterable[str] | None = None,
    spatial_key: str | None = None,
    output_dir: str | Path | None = None,
    filename_prefix: str | None = None,
    save_prefix: str | None = None,
    colors: dict[str, str] | None = None,
    palette: dict[str, str] | list[str] | tuple[str, ...] | str | None = None,
    fixed_window_um: float | None = None,
    center_method: str = "median",
    point_size: float = 4.0,
    underlay_size: float = 2.0,
    underlay_color: str = "#bdbdbd",
    underlay_alpha: float = 0.08,
    scale_bar: bool = True,
    scale_bar_um: float = 1000.0,
    scale_bar_label: str = "1 mm",
    figsize: tuple[float, float] = (10.0, 10.0),
    dpi: int = 600,
    legend_width: float = 2.2,
    combined: bool = False,
    save_individual: bool = True,
    max_cols: int = 3,
    show: bool = True,
    verbose: bool = True,
) -> dict[str, object]:
    """Create spatial cell type plots with a grey all-cell underlay.

    Set ``sample_col`` to use a shortened image label column such as ``ImageID``
    instead of the default QuPath ``Image`` column. By default, all selected
    samples are plotted in a shared centered window based on the largest sample
    x/y extent around each sample's median cell coordinate. Pass
    ``fixed_window_um`` to force a square fixed-size window.
    """

    plt, mtick, np, pd, Line2D, hsv_to_rgb, to_hex = _require_plotting()
    spatial_key = _resolve_spatial_key(adata, spatial_key)
    output_dir = _resolve_plot_dir(adata, output_dir) / "spatial_celltypes"
    output_dir.mkdir(parents=True, exist_ok=True)

    obs_plot = _prepare_obs(
        adata,
        category_col=category_col,
        sample_col=sample_col,
        subset_col=subset_col,
        subset_value=subset_value,
        samples=samples if samples is not None else images,
        categories=celltypes,
    )
    if obs_plot.empty:
        raise ValueError("No cells available for spatial plotting after filtering.")

    all_samples = sorted(adata.obs[sample_col].astype(str).unique().tolist())
    if samples is not None:
        selected_samples = list(samples)
    elif images is not None:
        selected_samples = list(images)
    else:
        selected_samples = all_samples
    selected_samples = [str(sample) for sample in selected_samples]
    if celltypes is not None:
        category_order = [str(celltype) for celltype in celltypes]
    else:
        category_order = sorted(obs_plot[category_col].astype(str).unique().tolist())
    _base_palette = _get_category_palette(adata, category_col)
    color_lookup = _color_map(
        category_order, hsv_to_rgb, to_hex, colors, palette if palette is not None else _base_palette
    )
    center_method = str(center_method).lower()
    if center_method not in {"median", "mean", "bbox"}:
        raise ValueError("center_method must be 'median', 'mean', or 'bbox'.")

    sample_bounds = {}
    for sample in selected_samples:
        sample_mask = adata.obs[sample_col].astype(str) == sample
        if not sample_mask.any():
            continue
        sample_positions = np.flatnonzero(sample_mask.to_numpy())
        sample_coords = adata.obsm[spatial_key][sample_positions, :]
        x_min = float(np.min(sample_coords[:, 0]))
        x_max = float(np.max(sample_coords[:, 0]))
        y_min = float(np.min(sample_coords[:, 1]))
        y_max = float(np.max(sample_coords[:, 1]))
        if center_method == "median":
            x_center = float(np.median(sample_coords[:, 0]))
            y_center = float(np.median(sample_coords[:, 1]))
        elif center_method == "mean":
            x_center = float(np.mean(sample_coords[:, 0]))
            y_center = float(np.mean(sample_coords[:, 1]))
        else:
            x_center = (x_min + x_max) / 2.0
            y_center = (y_min + y_max) / 2.0
        sample_bounds[sample] = {
            "x_min": x_min,
            "x_max": x_max,
            "y_min": y_min,
            "y_max": y_max,
            "x_center": x_center,
            "y_center": y_center,
            "width": x_max - x_min,
            "height": y_max - y_min,
            "x_radius": max(abs(x_min - x_center), abs(x_max - x_center)),
            "y_radius": max(abs(y_min - y_center), abs(y_max - y_center)),
        }

    if not sample_bounds:
        raise ValueError(f"No samples found for plotting using sample_col={sample_col!r}.")

    if fixed_window_um is not None:
        shared_width = float(fixed_window_um)
        shared_height = float(fixed_window_um)
    else:
        shared_width = 2.0 * max(bounds["x_radius"] for bounds in sample_bounds.values())
        shared_height = 2.0 * max(bounds["y_radius"] for bounds in sample_bounds.values())
        shared_width = max(shared_width, float(scale_bar_um) * 1.2, 1.0)
        shared_height = max(shared_height, 1.0)
        shared_width *= 1.04
        shared_height *= 1.04

    def _legend_handles():
        return [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                color=color_lookup[category],
                markersize=8,
                label=category,
            )
            for category in category_order
        ]

    def _draw_sample(ax, image: str) -> bool:
        image_mask = adata.obs[sample_col].astype(str) == image
        if not image_mask.any():
            return False
        image_positions = np.flatnonzero(image_mask.to_numpy())
        sub_obs = adata.obs.iloc[image_positions]
        coords = adata.obsm[spatial_key][image_positions, :]

        plot_mask = sub_obs.index.isin(obs_plot.index)
        if not plot_mask.any():
            return False
        plot_obs = sub_obs.loc[plot_mask]

        bounds = sample_bounds[image]
        x_half = shared_width / 2.0
        y_half = shared_height / 2.0
        x_lim = (-x_half, x_half)
        y_lim = (-y_half, y_half)
        centered_coords = coords.copy()
        centered_coords[:, 0] = centered_coords[:, 0] - bounds["x_center"]
        centered_coords[:, 1] = centered_coords[:, 1] - bounds["y_center"]
        plot_coords = centered_coords[plot_mask, :]

        ax.scatter(
            centered_coords[:, 0],
            centered_coords[:, 1],
            s=underlay_size,
            c=underlay_color,
            alpha=underlay_alpha,
            linewidths=0,
            marker="o",
            zorder=1,
        )

        values = plot_obs[category_col].astype(str).to_numpy()
        for category in category_order:
            mask = values == category
            if not mask.any():
                continue
            ax.scatter(
                plot_coords[mask, 0],
                plot_coords[mask, 1],
                s=point_size,
                c=[color_lookup[category]],
                alpha=1.0,
                linewidths=0,
                marker="o",
                zorder=3,
            )

        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        ax.set_aspect("equal")
        ax.set_anchor("C")
        ax.axis("off")
        if scale_bar:
            _add_scale_bar(ax, x_lim, y_lim, length_um=scale_bar_um, label=scale_bar_label)
        title = f"{subset_value} | {image}" if subset_value else image
        ax.set_title(title, fontsize=10, pad=10)
        return True

    fig_paths: list[Path] = []
    handles = _legend_handles()

    if save_individual:
        for image in selected_samples:
            fig = plt.figure(figsize=(figsize[0] + legend_width, figsize[1]), constrained_layout=False)
            gs = fig.add_gridspec(
                1,
                2,
                left=0.02,
                right=0.98,
                bottom=0.04,
                top=0.94,
                width_ratios=[figsize[0], legend_width],
                wspace=0.02,
            )
            ax = fig.add_subplot(gs[0, 0])
            legend_ax = fig.add_subplot(gs[0, 1])
            legend_ax.axis("off")
            plotted = _draw_sample(ax, image)
            if not plotted:
                plt.close(fig)
                continue
            if handles:
                legend_ax.legend(
                    handles=handles,
                    loc="center left",
                    frameon=False,
                    fontsize=8,
                    borderaxespad=0.0,
                )

            prefix_parts = ["spatial"]
            if subset_value:
                prefix_parts.append(_safe_name(subset_value))
            prefix_parts.append(_safe_name(image))
            if save_prefix:
                prefix = f"{save_prefix}_{_safe_name(image)}"
            else:
                prefix = filename_prefix or "_".join(prefix_parts)
            fig_path = output_dir / f"{prefix}.png"
            pdf_path = output_dir / f"{prefix}.pdf"
            fig.savefig(fig_path, dpi=dpi)
            fig.savefig(pdf_path, dpi=dpi)
            if show:
                plt.show()
            plt.close(fig)
            fig_paths.extend([fig_path, pdf_path])

            if verbose:
                print(f"Saved spatial plot: {fig_path}")

    combined_paths: list[Path] = []
    if combined:
        plotted_samples = [sample for sample in selected_samples if sample in sample_bounds]
        n_plots = len(plotted_samples)
        if n_plots:
            n_cols = max(1, min(int(max_cols), n_plots))
            n_rows = int(np.ceil(n_plots / n_cols))
            fig = plt.figure(
                figsize=((figsize[0] / 2) * n_cols + legend_width, (figsize[1] / 2) * n_rows),
                constrained_layout=False,
            )
            gs = fig.add_gridspec(
                n_rows,
                n_cols + 1,
                left=0.02,
                right=0.98,
                bottom=0.04,
                top=0.94,
                width_ratios=[*[1] * n_cols, legend_width / 5],
                wspace=0.02,
                hspace=0.08,
            )
            for index, sample in enumerate(plotted_samples):
                row = index // n_cols
                col = index % n_cols
                ax = fig.add_subplot(gs[row, col])
                if not _draw_sample(ax, sample):
                    ax.axis("off")
            for index in range(n_plots, n_rows * n_cols):
                row = index // n_cols
                col = index % n_cols
                fig.add_subplot(gs[row, col]).axis("off")
            legend_ax = fig.add_subplot(gs[:, -1])
            legend_ax.axis("off")
            if handles:
                legend_ax.legend(
                    handles=handles,
                    loc="center left",
                    frameon=False,
                    fontsize=8,
                    borderaxespad=0.0,
                )
            prefix_parts = ["spatial_combined"]
            if subset_value:
                prefix_parts.append(_safe_name(subset_value))
            prefix = save_prefix or filename_prefix or "_".join(prefix_parts)
            fig_path = output_dir / f"{prefix}.png"
            pdf_path = output_dir / f"{prefix}.pdf"
            fig.savefig(fig_path, dpi=dpi)
            fig.savefig(pdf_path, dpi=dpi)
            if show:
                plt.show()
            plt.close(fig)
            combined_paths.extend([fig_path, pdf_path])
            fig_paths.extend(combined_paths)
            if verbose:
                print(f"Saved combined spatial plot: {fig_path}")

    return {
        "figures": fig_paths,
        "combined_figures": combined_paths,
        "shared_extent_um": {"width": shared_width, "height": shared_height},
        "center_method": center_method,
        "sample_bounds": sample_bounds,
    }


# ── Heatmap infrastructure ───────────────────────────────────────────────────
#
# Tile size:  0.25 cm × 0.25 cm per heatmap cell
# Journal widths (width=):
#   "single" → 90 mm  (one journal column)
#   "double" → 180 mm (two journal columns / full page width)
#   "auto"   → tile-based (explicit opt-in only)
#
# Saved formats: PDF (vector), SVG (vector), TIFF (dpi-raster, default 300)
# ─────────────────────────────────────────────────────────────────────────────

try:
    import cmcrameri  # noqa: F401
    _CRAMERI_AVAILABLE = True
    _DEFAULT_SEQ_CMAP  = "cmc.batlow"
except ImportError:
    _CRAMERI_AVAILABLE = False
    _DEFAULT_SEQ_CMAP  = "cividis"

# Short-name aliases for Crameri colormaps with fallbacks for when cmcrameri
# is not installed. Pass these names as cmap= to any heatmap function.
# e.g. cmap="roma", cmap="batlow", cmap="vik"
_CRAMERI_ALIASES: "dict[str, tuple[str, str]]" = {
    # name      → (cmc name,        fallback)
    "batlow":    ("cmc.batlow",     "cividis"),
    "batlowS":   ("cmc.batlowS",    "tab20"),
    "roma":      ("cmc.roma",       "coolwarm"),
    "vik":       ("cmc.vik",        "coolwarm"),
    "berlin":    ("cmc.berlin",     "RdBu_r"),
    "lisbon":    ("cmc.lisbon",     "RdBu"),
    "cork":      ("cmc.cork",       "PiYG"),
    "broc":      ("cmc.broc",       "PRGn"),
    "oleron":    ("cmc.oleron",     "terrain"),
    "nuuk":      ("cmc.nuuk",       "Blues"),
    "lapaz":     ("cmc.lapaz",      "plasma"),
    "tokyo":     ("cmc.tokyo",      "magma"),
}


def _resolve_cmap(name: str) -> str:
    """Resolve a colormap name, expanding Crameri short aliases.

    If *name* is a key in ``_CRAMERI_ALIASES``, returns the ``cmc.*`` name
    when ``cmcrameri`` is installed, or the fallback matplotlib name otherwise.
    Unknown names are returned as-is (passed directly to matplotlib).
    """
    if name in _CRAMERI_ALIASES:
        cmc_name, fallback = _CRAMERI_ALIASES[name]
        return cmc_name if _CRAMERI_AVAILABLE else fallback
    return name

_MM_TO_IN      = 1.0 / 25.4
_TILE_MM       = 4.0                         # 0.25 cm per cell
_TILE_IN       = _TILE_MM * _MM_TO_IN        # ≈ 0.0984 in

_SINGLE_COL_IN = 90.0  * _MM_TO_IN          # 90 mm  ≈ 3.54 in
_DOUBLE_COL_IN = 180.0 * _MM_TO_IN          # 180 mm ≈ 7.09 in

# Internal layout constants (inches)
_STRIP_W   = 0.15   # width/height of one annotation colour strip
_STRIP_GAP = 0.03   # gap between consecutive strips
_CBAR_W    = 0.18   # colorbar width
_CBAR_PAD  = 0.10   # gap: heatmap right edge → colorbar
_TITLE_H   = 0.32   # vertical space above colour strips for title text
_BOT_PAD   = 0.08   # bottom figure margin (col-label text extends below via bbox_inches)


def _glasbey_colors(n: int) -> list:
    """Return *n* maximally-distinct colours.

    Resolution order: ``glasbey`` package → ``colorcet`` glasbey palette →
    matplotlib ``tab20`` / ``hsv`` fallback.
    """
    if n <= 0:
        return []
    try:
        import glasbey as _gb
        return list(_gb.create_palette(palette_size=max(n, 2)))[:n]
    except ImportError:
        pass
    try:
        import colorcet as _cc
        pal = _cc.palette["glasbey_bw_minc_20"]
        return [pal[i % len(pal)] for i in range(n)]
    except (ImportError, KeyError):
        pass
    import matplotlib.pyplot as _plt
    base = _plt.get_cmap("tab20" if n <= 20 else "hsv")
    return [base(i / max(n, 1)) for i in range(n)]


def _cluster_order(matrix) -> "tuple[list[int], list[int]]":
    """Average-linkage hierarchical clustering.

    Returns ``(row_order, col_order)`` leaf-index lists for *matrix*.
    Both are computed on the *original* matrix (independent permutations).
    Requires scipy.
    """
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist
    except ImportError as exc:
        raise ImportError(
            "scipy is required for heatmap clustering. "
            "Install it with:  pip install scipy"
        ) from exc
    import numpy as _np

    def _leaves(m: "_np.ndarray") -> "list[int]":
        if m.shape[0] < 2:
            return list(range(m.shape[0]))
        d = pdist(_np.nan_to_num(m, nan=0.0), metric="euclidean")
        return list(leaves_list(linkage(d, method="average")))

    return _leaves(matrix), _leaves(matrix.T)


def _cat_strip_colors(labels: list, palette: "dict | None" = None) -> "tuple[dict, list]":
    """Build ``{label: colour}`` mapping using the glasbey palette.

    If *palette* is provided (a pre-built ``{label: hex}`` dict), colours from
    it are used preferentially so strips stay consistent with other plots.
    Returns ``(color_map, ordered_categories)``.
    """
    import matplotlib.colors as _mcolors

    cats = list(dict.fromkeys(str(l) for l in labels))   # unique, insertion-ordered
    if palette:
        out: dict[str, str] = {}
        missing = []
        for c in cats:
            if c in palette:
                col = palette[c]
                out[c] = col if isinstance(col, str) else _mcolors.to_hex(col)
            else:
                missing.append(c)
        if missing:
            extra = _glasbey_colors(len(missing))
            for c, col in zip(missing, extra):
                out[c] = col if isinstance(col, str) else _mcolors.to_hex(col)
        return out, cats
    colors = _glasbey_colors(len(cats))
    return {c: colors[i] for i, c in enumerate(cats)}, cats


def _build_heatmap_figure(
    matrix,
    row_labels: list,
    col_labels: list,
    *,
    # colour scale
    cmap: str,
    vmin: "float | None",
    vmax: "float | None",
    center: "float | None",
    # cell annotations
    annotate: bool,
    fmt: str,
    # text
    title: str,
    cbar_label: str,
    # clustering (no dendrograms)
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    # annotation strips
    # each strip: (label_list_in_original_order, strip_name_or_None)
    row_strips: "list | None" = None,   # drawn left of heatmap
    col_strips: "list | None" = None,   # drawn above heatmap
    row_strip_palette: "dict | None" = None,  # pre-built {label: colour} for row strips
    # layout
    width: str = "auto",                # "auto" | "single" | "double"
) -> "tuple":
    """Compose a publication-ready heatmap figure.

    Returns ``(matplotlib.figure.Figure, ordered_matrix_DataFrame)``.
    Row/column order reflects clustering (if enabled).
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    import pandas as pd

    row_strips = list(row_strips or [])
    col_strips = list(col_strips or [])
    n_r, n_c   = matrix.shape

    # ── 1. Hierarchical clustering (reorder rows and cols, no dendrograms) ──
    row_ord = list(range(n_r))
    col_ord = list(range(n_c))
    if (cluster_rows and n_r > 1) or (cluster_cols and n_c > 1):
        r_all, c_all = _cluster_order(matrix)
        if cluster_rows and n_r > 1:
            row_ord = r_all
        if cluster_cols and n_c > 1:
            col_ord = c_all

    matrix     = matrix[np.ix_(row_ord, col_ord)]
    row_labels = [row_labels[i] for i in row_ord]
    col_labels = [col_labels[i] for i in col_ord]
    row_strips = [([s[0][i] for i in row_ord], s[1]) for s in row_strips]
    col_strips = [([s[0][i] for i in col_ord], s[1]) for s in col_strips]

    # ── 2. Colour scale ─────────────────────────────────────────────────────
    if center is not None:
        vmin = vmin if vmin is not None else -3.0
        vmax = vmax if vmax is not None else  3.0

    # ── 3. Figure geometry (all in inches) ──────────────────────────────────
    n_rs = len(row_strips)
    n_cs = len(col_strips)
    rstrip_w = n_rs * (_STRIP_W + _STRIP_GAP)
    cstrip_h = n_cs * (_STRIP_W + _STRIP_GAP)

    core_h = n_r * _TILE_IN
    core_w = n_c * _TILE_IN

    # Estimate width needed for row label text (right side of heatmap).
    # At 7 pt sans-serif each character is ≈ 4.2 pt wide; add 4 pt tick pad.
    _max_row_chars = max((len(str(l)) for l in row_labels), default=4)
    _row_lbl_w = max((_max_row_chars * 4.2 + 4) / 72, 0.35)  # inches

    # Height: bottom pad + core + col strips above + title space
    # (rotated col-label text extends below BOT_PAD; captured by bbox_inches="tight")
    fig_h = _BOT_PAD + core_h + cstrip_h + _TITLE_H

    _right_margin = _row_lbl_w + _CBAR_PAD + _CBAR_W   # heatmap → right edge

    # Tile width is always fixed at _TILE_IN (3.5 mm); figure width follows from
    # the tile count. width="single"/"double" set a minimum figure width — if
    # the tile-based width already exceeds it, the tile size wins.
    core_w_use = core_w   # always tile-based
    tile_fig_w = rstrip_w + core_w_use + _right_margin

    if width == "single":
        fig_w = max(tile_fig_w, _SINGLE_COL_IN)
    elif width == "double":
        fig_w = max(tile_fig_w, _DOUBLE_COL_IN)
    else:  # "auto"
        fig_w = tile_fig_w

    # Axes coordinates in figure-fraction units (origin = bottom-left)
    def _fx(v: float) -> float: return v / fig_w   # noqa: E731
    def _fy(v: float) -> float: return v / fig_h   # noqa: E731

    y0  = _fy(_BOT_PAD)    # bottom of heatmap core
    h_c = _fy(core_h)      # heatmap core height (fraction)
    x0  = _fx(rstrip_w)    # left edge of heatmap core
    w_c = _fx(core_w_use)  # heatmap core width (fraction)

    fig = plt.figure(figsize=(fig_w, fig_h))

    # ── 4. Main heatmap ─────────────────────────────────────────────────────
    ax_h = fig.add_axes([x0, y0, w_c, h_c])
    im = ax_h.imshow(
        matrix, aspect="auto", cmap=cmap,
        vmin=vmin, vmax=vmax, interpolation="nearest",
    )
    # Col labels (bottom)
    ax_h.set_xticks(range(n_c))
    ax_h.set_xticklabels(col_labels, rotation=90, ha="center", fontsize=7)
    ax_h.xaxis.tick_bottom()
    ax_h.tick_params(axis="x", length=0, pad=2)
    # Row labels (right)
    ax_h.set_yticks(range(n_r))
    ax_h.set_yticklabels(row_labels, fontsize=7)
    ax_h.yaxis.tick_right()
    ax_h.tick_params(axis="y", length=0, pad=2)
    for sp in ax_h.spines.values():
        sp.set_visible(False)

    # Cell annotations with luminance-based text colour
    if annotate:
        _vmin_ann = vmin if vmin is not None else float(np.nanmin(matrix))
        _vmax_ann = vmax if vmax is not None else float(np.nanmax(matrix))
        _norm_ann = plt.Normalize(vmin=_vmin_ann, vmax=_vmax_ann)
        _cmap_obj = plt.get_cmap(cmap)
        for r in range(n_r):
            for c in range(n_c):
                v = matrix[r, c]
                if np.isnan(v):
                    continue
                bg_rgba = _cmap_obj(_norm_ann(v))
                # Perceived luminance (ITU-R BT.601)
                lum = 0.299 * bg_rgba[0] + 0.587 * bg_rgba[1] + 0.114 * bg_rgba[2]
                tc  = "white" if lum < 0.50 else "black"
                ax_h.text(c, r, f"{v:{fmt}}", ha="center", va="center",
                           fontsize=5.5, color=tc, clip_on=True)

    # ── 5. Colorbar ─────────────────────────────────────────────────────────
    x_cb = _fx(rstrip_w + core_w_use + _row_lbl_w + _CBAR_PAD)
    ax_cb = fig.add_axes([x_cb, y0, _fx(_CBAR_W), h_c])
    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("")
    cb.ax.set_xlabel(cbar_label, fontsize=7, labelpad=4, ha="center")
    cb.ax.tick_params(labelsize=5, length=2, width=0.5)
    for sp in cb.ax.spines.values():
        sp.set_visible(False)

    # ── 6. Title (above col strips) ─────────────────────────────────────────
    t_y = y0 + h_c + _fy(cstrip_h) + _fy(_TITLE_H * 0.15)
    fig.text(x0, t_y, title, fontsize=8, va="bottom", ha="left")

    # ── 7. Row strips (left of heatmap; row 0 at top to match imshow) ───────
    for si, (strip_lbls, _strip_name) in enumerate(row_strips):
        color_map, _ = _cat_strip_colors(strip_lbls, palette=row_strip_palette)
        sx = _fx(si * (_STRIP_W + _STRIP_GAP))
        ax_s = fig.add_axes([sx, y0, _fx(_STRIP_W), h_c])
        for ri, lbl in enumerate(strip_lbls):
            ax_s.add_patch(mpatches.Rectangle(
                (0, ri - 0.5), 1, 1,
                facecolor=color_map[str(lbl)], linewidth=0,
            ))
        ax_s.set_xlim(0, 1)
        ax_s.set_ylim(n_r - 0.5, -0.5)   # inverted: row 0 at top
        ax_s.axis("off")

    # ── 8. Col strips (above heatmap) + legend ──────────────────────────────
    legend_patches = []
    for si, (strip_lbls, strip_name) in enumerate(col_strips):
        color_map, cats = _cat_strip_colors(strip_lbls)
        sy = y0 + h_c + _fy(si * (_STRIP_W + _STRIP_GAP))
        ax_s = fig.add_axes([x0, sy, w_c, _fy(_STRIP_W)])
        for ci, lbl in enumerate(strip_lbls):
            ax_s.add_patch(mpatches.Rectangle(
                (ci - 0.5, 0), 1, 1,
                facecolor=color_map[str(lbl)], linewidth=0,
            ))
        ax_s.set_xlim(-0.5, n_c - 0.5)
        ax_s.set_ylim(0, 1)
        ax_s.axis("off")
        lbl_prefix = f"{strip_name}: " if strip_name else ""
        for cat in cats:
            legend_patches.append(
                mpatches.Patch(facecolor=color_map[cat],
                               label=f"{lbl_prefix}{cat}", linewidth=0)
            )

    # Legend for col-strip categories (placed outside figure; captured by bbox_inches="tight")
    if legend_patches:
        fig.legend(
            handles=legend_patches,
            fontsize=6, frameon=True, framealpha=0.9,
            edgecolor="#cccccc", ncol=1,
            loc="lower right", bbox_to_anchor=(1.0, 0.0),
            borderpad=0.5, handlelength=1.0,
        )

    matrix_df = pd.DataFrame(matrix, index=row_labels, columns=col_labels)
    return fig, matrix_df


# ── Public heatmap functions ─────────────────────────────────────────────────

def plot_marker_heatmap(
    adata,
    *,
    category_col: str = "celltype",
    markers: "list[str] | None" = None,
    values: str = "positivity",
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    width: str = "single",
    cmap: "str | None" = None,
    annotate: bool = False,
    row_strip: bool = False,
    dpi: int = 600,
    output_dir: "str | Path | None" = None,
    show: bool = True,
    verbose: bool = True,
) -> "dict[str, list]":
    """Plot a heatmap of marker expression per cell type (or any category).

    Rows are clustered by hierarchical clustering (no dendrogram shown).
    Columns (markers) are also clustered independently.
    A colour strip on the left uses the glasbey palette to colour each row category.

    Parameters
    ----------
    adata:
        AnnData object with ``adata.obs[category_col]`` present.
    category_col:
        Column in ``adata.obs`` to group cells by (default ``"celltype"``).
    markers:
        List of marker names to include. Defaults to all markers in
        ``adata.var_names``.
    values:
        ``"positivity"`` — mean fraction of cells positive (from ``_pos``
        columns in ``adata.obs``); colormap defaults to ``"cividis"``,
        scale 0–1.

        ``"intensity"``  — mean marker intensity from ``adata.X``,
        Z-score normalised per marker; colormap defaults to ``"coolwarm"``,
        centred at 0, clipped ±3.

        ``"both"``       — generate both plots in one call.
    cluster_rows:
        Reorder rows by hierarchical clustering (default True).
    cluster_cols:
        Reorder columns by hierarchical clustering (default True).
    width:
        Figure width preset.
        ``"single"`` — 90 mm, one journal column (default).
        ``"double"`` — 180 mm, two journal columns / full page width.
        ``"auto"``   — tile-based sizing (0.25 cm × 0.25 cm per cell); produces very small figures for typical panel sizes — use only when you need exact tile dimensions.
    cmap:
        Override the default colormap. Pass any matplotlib colormap name.
    annotate:
        Write values inside each cell (default True).
    dpi:
        Resolution for the TIFF output (default 300). PDF and SVG are
        vector formats and ignore this parameter.
    output_dir:
        Folder to save plots. Defaults to the QXYCell output folder.
    show:
        Display the figure inline (default True).
    verbose:
        Print save paths (default True).

    Returns
    -------
    dict mapping mode → ``[pdf_path, svg_path, tif_path, csv_path]``.
    """
    plt, mtick, np, pd, Line2D, hsv_to_rgb, to_hex = _require_plotting()
    out_dir = _resolve_plot_dir(adata, output_dir) / "marker_heatmap"
    out_dir.mkdir(parents=True, exist_ok=True)

    if category_col not in adata.obs.columns:
        raise ValueError(f"'{category_col}' not found in adata.obs.")

    # Pre-build (or retrieve) the global glasbey palette for this category column
    # so the row strip uses the same colours as plot_spatial / plot_stacked_bar.
    _cat_palette = _get_category_palette(adata, category_col) if row_strip else None

    marker_names = list(markers) if markers is not None else list(adata.var_names)
    modes = ["positivity", "intensity"] if values == "both" else [values]
    saved: dict[str, list] = {}

    for mode in modes:
        if mode == "positivity":
            pos_cols = [f"{m}_pos" for m in marker_names
                        if f"{m}_pos" in adata.obs.columns]
            if not pos_cols:
                raise ValueError(
                    "No '_pos' columns found. Run qxy.run() first."
                )
            display_markers = [c.removesuffix("_pos") for c in pos_cols]
            matrix_df  = adata.obs.groupby(category_col, observed=True)[pos_cols].mean()
            matrix_df.columns = display_markers
            matrix     = matrix_df.values.astype(float)
            row_labels = list(matrix_df.index)
            col_labels = display_markers
            cbar_label = "f"
            _cmap      = _resolve_cmap(cmap) if cmap is not None else _DEFAULT_SEQ_CMAP
            _vmin, _vmax, _center, _fmt = 0.0, 1.0, None, ".2f"
            title      = f"Marker positivity by {category_col}"

        elif mode == "intensity":
            X = adata.X
            # Handle sparse matrices without hard-importing scipy at module level
            if hasattr(X, "toarray"):
                X = X.toarray()
            X = np.asarray(X, dtype=float)
            valid = [m for m in marker_names if m in adata.var_names]
            if not valid:
                raise ValueError(
                    "None of the specified markers found in adata.var_names."
                )
            midx  = [list(adata.var_names).index(m) for m in valid]
            X_sub = X[:, midx]
            obs_df = adata.obs[[category_col]].copy()
            obs_df["_idx"] = range(len(obs_df))
            groups = obs_df.groupby(category_col, observed=True)["_idx"].apply(list)
            mat = np.zeros((len(groups), len(valid)))
            for r, (_, idxs) in enumerate(groups.items()):
                mat[r] = X_sub[idxs].mean(axis=0)
            col_std = mat.std(axis=0)
            col_std[col_std == 0] = 1
            matrix     = (mat - mat.mean(axis=0)) / col_std
            row_labels = list(groups.index)
            col_labels = valid
            matrix_df  = pd.DataFrame(matrix, index=row_labels, columns=col_labels)
            cbar_label = "z (↓)"
            _cmap      = _resolve_cmap(cmap) if cmap is not None else "coolwarm"
            _vmin, _vmax, _center, _fmt = -3.0, 3.0, 0.0, ".2f"
            title      = f"Marker intensity by {category_col} (Z-score)"

        else:
            raise ValueError("values must be 'positivity', 'intensity', or 'both'.")

        row_strips = [(list(row_labels), None)] if row_strip else []

        fig, out_df = _build_heatmap_figure(
            matrix, row_labels, col_labels,
            cmap=_cmap, vmin=_vmin, vmax=_vmax, center=_center,
            annotate=annotate, fmt=_fmt,
            title=title, cbar_label=cbar_label,
            cluster_rows=cluster_rows, cluster_cols=cluster_cols,
            row_strips=row_strips, col_strips=[],
            row_strip_palette=_cat_palette,
            width=width,
        )

        prefix   = f"marker_heatmap_{mode}"
        pdf_path = out_dir / f"{prefix}.pdf"
        svg_path = out_dir / f"{prefix}.svg"
        tif_path = out_dir / f"{prefix}.tif"
        csv_path = out_dir / f"{prefix}.csv"
        fig.savefig(pdf_path, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight", format="svg")
        fig.savefig(tif_path, dpi=dpi, bbox_inches="tight", format="tiff")
        out_df.to_csv(csv_path)
        if show:
            plt.show()
        plt.close(fig)
        saved[mode] = [pdf_path, svg_path, tif_path, csv_path]
        if verbose:
            print(f"Saved marker heatmap ({mode}): {pdf_path.parent}/")

    return saved


def plot_cn_heatmap(
    adata,
    *,
    cn_col: str = "cn",
    sample_col: str = "Image",
    condition_col: "str | None" = None,
    normalize: str = "sample",
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    width: str = "single",
    cmap: "str | None" = None,
    annotate: bool = False,
    row_strip: bool = False,
    dpi: int = 600,
    output_dir: "str | Path | None" = None,
    show: bool = True,
    verbose: bool = True,
) -> "dict[str, list]":
    """Plot a heatmap of CN abundance across samples.

    Rows (CNs) and columns (samples) are clustered by hierarchical clustering
    (no dendrogram shown). A colour strip on the left uses the glasbey palette
    to colour each CN row. If ``condition_col`` is provided, a colour strip
    above the heatmap groups samples by their condition.

    Parameters
    ----------
    adata:
        AnnData object with ``adata.obs[cn_col]`` and ``adata.obs[sample_col]``.
    cn_col:
        Column in ``adata.obs`` containing CN labels (default ``"cn"``).
    sample_col:
        Column in ``adata.obs`` containing sample labels (default ``"Image"``).
    condition_col:
        Optional column in ``adata.obs`` mapping samples to experimental
        conditions. When provided, a coloured annotation strip is drawn above
        the heatmap grouping samples by condition, with a legend.
    normalize:
        ``"sample"`` — each column (sample) sums to 1; shows CN composition
        within each sample.

        ``"cn"``     — each row (CN) sums to 1; shows sample origin of each CN.

        ``"both"``   — generate both normalisations in one call.
    cluster_rows:
        Reorder CN rows by hierarchical clustering (default True).
    cluster_cols:
        Reorder sample columns by hierarchical clustering (default True).
    width:
        Figure width preset.
        ``"single"`` — 90 mm, one journal column (default).
        ``"double"`` — 180 mm, two journal columns / full page width.
        ``"auto"``   — tile-based sizing (0.25 cm × 0.25 cm per cell); produces very small figures for typical panel sizes — use only when you need exact tile dimensions.
    cmap:
        Colormap name. Defaults to ``"cividis"``.
    annotate:
        Write values inside each cell (default True).
    dpi:
        Resolution for the TIFF output (default 300). PDF and SVG ignore this.
    output_dir:
        Folder to save plots. Defaults to the QXYCell output folder.
    show:
        Display the figure inline (default True).
    verbose:
        Print save paths (default True).

    Returns
    -------
    dict mapping normalize mode → ``[pdf_path, svg_path, tif_path, csv_path]``.
    """
    plt, mtick, np, pd, Line2D, hsv_to_rgb, to_hex = _require_plotting()
    out_dir = _resolve_plot_dir(adata, output_dir) / "cn_heatmap"
    out_dir.mkdir(parents=True, exist_ok=True)

    if cn_col not in adata.obs.columns:
        raise ValueError(
            f"'{cn_col}' not found in adata.obs. Run qxy.cn_kmeans() first."
        )
    if sample_col not in adata.obs.columns:
        raise ValueError(f"'{sample_col}' not found in adata.obs.")
    if condition_col is not None and condition_col not in adata.obs.columns:
        raise ValueError(f"'{condition_col}' not found in adata.obs.")

    _cmap = _resolve_cmap(cmap) if cmap is not None else _DEFAULT_SEQ_CMAP
    modes = ["sample", "cn"] if normalize == "both" else [normalize]
    saved: dict[str, list] = {}

    # Pre-build (or retrieve) the shared CN palette so the row strip uses the
    # same colours as plot_stacked_bar(category_col="cn") and plot_spatial(category_col="cn").
    _cn_palette = _get_category_palette(adata, cn_col) if row_strip else None

    for mode in modes:
        norm_arg = "columns" if mode == "sample" else "index"
        ct = pd.crosstab(
            adata.obs[cn_col].astype(str),
            adata.obs[sample_col].astype(str),
            normalize=norm_arg,
        )
        matrix     = ct.values.astype(float)
        row_labels = list(ct.index)    # CNs
        col_labels = list(ct.columns)  # samples

        if mode == "sample":
            title      = f"CN fraction per {sample_col} (columns sum to 1)"
            cbar_label = "f (↓)"
        else:
            title      = f"CN composition by {sample_col} (rows sum to 1)"
            cbar_label = "f (→)"

        row_strips = [(list(row_labels), None)] if row_strip else []

        # Top strip: sample → condition (if condition_col supplied)
        col_strips: list = []
        if condition_col is not None:
            cond_lookup = (
                adata.obs[[sample_col, condition_col]]
                .drop_duplicates(subset=[sample_col])
                .set_index(sample_col)[condition_col]
            )
            # Aligned to ct.columns (original order; clustering reorders inside _build)
            cond_labels = [str(cond_lookup.get(s, "Unknown")) for s in col_labels]
            col_strips  = [(cond_labels, condition_col)]

        fig, out_df = _build_heatmap_figure(
            matrix, row_labels, col_labels,
            cmap=_cmap, vmin=0.0, vmax=None, center=None,
            annotate=annotate, fmt=".2f",
            title=title, cbar_label=cbar_label,
            cluster_rows=cluster_rows, cluster_cols=cluster_cols,
            row_strips=row_strips, col_strips=col_strips,
            row_strip_palette=_cn_palette,
            width=width,
        )

        prefix   = f"cn_heatmap_by_{mode}"
        pdf_path = out_dir / f"{prefix}.pdf"
        svg_path = out_dir / f"{prefix}.svg"
        tif_path = out_dir / f"{prefix}.tif"
        csv_path = out_dir / f"{prefix}.csv"
        fig.savefig(pdf_path, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight", format="svg")
        fig.savefig(tif_path, dpi=dpi, bbox_inches="tight", format="tiff")
        out_df.to_csv(csv_path)
        if show:
            plt.show()
        plt.close(fig)
        saved[mode] = [pdf_path, svg_path, tif_path, csv_path]
        if verbose:
            print(f"Saved CN heatmap (normalize='{mode}'): {pdf_path.parent}/")

    return saved
