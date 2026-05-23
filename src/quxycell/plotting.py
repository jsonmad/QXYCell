"""Generic plotting helpers for QUXYCell AnnData objects."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Iterable


DEFAULT_PLOT_DIR = Path("outputs") / "qxy_run" / "plots"


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
            "QUXYCell plotting requires plot dependencies. Install with "
            "`pip install -e '.[plot]'` or `pip install quxycell[plot]`."
        ) from exc
    return plt, mtick, np, pd, Line2D, hsv_to_rgb, to_hex


def _safe_name(value: object) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(value)).strip("_")
    return safe or "plot"


def _fallback_color(label: object, hsv_to_rgb, to_hex) -> str:
    digest = hashlib.sha1(str(label).encode("utf-8")).hexdigest()
    hue = (int(digest[:8], 16) % 360) / 360.0
    return to_hex(hsv_to_rgb((hue, 0.72, 0.90)))


def _color_map(labels: Iterable[object], hsv_to_rgb, to_hex, colors: dict[str, str] | None = None):
    colors = colors or {}
    out = {}
    for label in labels:
        key = str(label)
        out[key] = colors.get(key, _fallback_color(key, hsv_to_rgb, to_hex))
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
    output_dir: str | Path = DEFAULT_PLOT_DIR,
    filename_prefix: str | None = None,
    colors: dict[str, str] | None = None,
    denominator: str = "all_cells",
    figsize: tuple[float, float] = (8.0, 4.0),
    dpi: int = 300,
    legend_width: float = 3.0,
    show: bool = True,
    verbose: bool = True,
) -> dict[str, Path]:
    """Create a generic stacked bar plot from AnnData cell annotations.

    Frequencies are computed per sample column, with optional aggregation by
    ``group_col``. Set ``sample_col`` to use a shortened image label column such
    as ``ImageID`` instead of the default QuPath ``Image`` column. By default,
    each category count is divided by all cells in the sample.
    """

    plt, mtick, np, pd, Line2D, hsv_to_rgb, to_hex = _require_plotting()
    output_dir = Path(output_dir).expanduser().resolve() / "stacked_bar"
    output_dir.mkdir(parents=True, exist_ok=True)
    table_dir = output_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    obs_plot = _prepare_obs(
        adata,
        category_col=category_col,
        sample_col=sample_col,
        group_col=group_col,
        subset_col=subset_col,
        subset_value=subset_value,
    )
    if obs_plot.empty:
        raise ValueError("No cells available for stacked bar plotting after filtering.")

    obs_all = adata.obs.copy()
    obs_all[sample_col] = obs_all[sample_col].astype(str).str.strip()
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
    prefix = filename_prefix or "_".join(safe_parts)

    color_lookup = _color_map(category_order, hsv_to_rgb, to_hex, colors)
    stack_colors = [color_lookup[str(label)] for label in category_order]

    fig = plt.figure(figsize=(figsize[0] + legend_width, figsize[1]), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[figsize[0], legend_width], wspace=0.08)
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
    ax.set_title(subset_value if subset_value else "Cell type frequency")
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.get_legend()
    if leg is not None:
        leg.remove()
    legend_ax.legend(
        handles,
        labels,
        loc="center left",
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

    table_path = table_dir / f"{prefix}_frequency.csv"
    freq_plot.to_csv(table_path)
    per_image_path = table_dir / f"{prefix}_per_image_frequency.csv"
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


def plot_spatial_celltypes(
    adata,
    *,
    category_col: str = "celltype",
    sample_col: str = "Image",
    subset_col: str | None = None,
    subset_value: str | None = None,
    images: Iterable[str] | None = None,
    spatial_key: str | None = None,
    output_dir: str | Path = DEFAULT_PLOT_DIR,
    filename_prefix: str | None = None,
    colors: dict[str, str] | None = None,
    fixed_window_um: float | None = None,
    center_method: str = "median",
    point_size: float = 4.0,
    underlay_size: float = 2.0,
    underlay_color: str = "#bdbdbd",
    underlay_alpha: float = 0.08,
    scale_bar_um: float = 1000.0,
    scale_bar_label: str = "1 mm",
    dpi: int = 300,
    legend_width: float = 2.2,
    show: bool = True,
    verbose: bool = True,
) -> dict[str, list[Path]]:
    """Create spatial cell type plots with a grey all-cell underlay.

    Set ``sample_col`` to use a shortened image label column such as ``ImageID``
    instead of the default QuPath ``Image`` column. By default, all selected
    samples are plotted in a shared centered window based on the largest sample
    x/y extent around each sample's median cell coordinate. Pass
    ``fixed_window_um`` to force a square fixed-size window.
    """

    plt, mtick, np, pd, Line2D, hsv_to_rgb, to_hex = _require_plotting()
    spatial_key = _resolve_spatial_key(adata, spatial_key)
    output_dir = Path(output_dir).expanduser().resolve() / "spatial_celltypes"
    output_dir.mkdir(parents=True, exist_ok=True)

    obs_plot = _prepare_obs(
        adata,
        category_col=category_col,
        sample_col=sample_col,
        subset_col=subset_col,
        subset_value=subset_value,
    )
    if obs_plot.empty:
        raise ValueError("No cells available for spatial plotting after filtering.")

    all_samples = sorted(adata.obs[sample_col].astype(str).unique().tolist())
    selected_samples = list(images) if images is not None else all_samples
    selected_samples = [str(sample) for sample in selected_samples]
    category_order = sorted(obs_plot[category_col].astype(str).unique().tolist())
    color_lookup = _color_map(category_order, hsv_to_rgb, to_hex, colors)
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

    fig_paths: list[Path] = []
    for image in selected_samples:
        image_mask = adata.obs[sample_col].astype(str) == image
        if not image_mask.any():
            continue
        image_positions = np.flatnonzero(image_mask.to_numpy())
        sub_obs = adata.obs.iloc[image_positions]
        coords = adata.obsm[spatial_key][image_positions, :]

        plot_mask = sub_obs.index.isin(obs_plot.index)
        if not plot_mask.any():
            continue
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

        fig = plt.figure(figsize=(10 + legend_width, 10), constrained_layout=False)
        gs = fig.add_gridspec(
            1,
            2,
            left=0.02,
            right=0.98,
            bottom=0.04,
            top=0.94,
            width_ratios=[10, legend_width],
            wspace=0.02,
        )
        ax = fig.add_subplot(gs[0, 0])
        legend_ax = fig.add_subplot(gs[0, 1])
        legend_ax.axis("off")
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
        _add_scale_bar(ax, x_lim, y_lim, length_um=scale_bar_um, label=scale_bar_label)
        title = f"{subset_value} | {image}" if subset_value else image
        ax.set_title(title, fontsize=10)

        handles = [
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

    return {
        "figures": fig_paths,
        "shared_extent_um": {"width": shared_width, "height": shared_height},
        "center_method": center_method,
        "sample_bounds": sample_bounds,
    }
