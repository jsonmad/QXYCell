"""QC summary reports for QXYCell AnnData objects."""

from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Any

from qxycell.paths import resolve_output_dir


def _write_table(table, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)
    return path


def qc(
    adata,
    *,
    sample_col: str = "Image",
    celltype_col: str = "celltype",
    annotation_prefix: str = "annotation__",
    marker_suffix: str = "_pos",
    output_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a compact QC report for a QXYCell AnnData object."""

    import pandas as pd

    out_dir = resolve_output_dir(output_dir, adata=adata) / "qc"
    out_dir.mkdir(parents=True, exist_ok=True)

    obs = adata.obs
    paths: dict[str, Path] = {}
    tables: dict[str, Any] = {}

    sample_values = obs[sample_col].astype(str) if sample_col in obs.columns else pd.Series([], dtype=str)
    annotation_cols = [column for column in obs.columns if str(column).startswith(annotation_prefix)]
    marker_cols = [column for column in obs.columns if str(column).endswith(marker_suffix)]
    ignore_cols = [
        column
        for column in annotation_cols
        if "ignore" in str(column).lower()
    ]

    overview = pd.DataFrame(
        [
            {"metric": "n_cells", "value": int(adata.n_obs)},
            {"metric": "n_markers", "value": int(adata.n_vars)},
            {"metric": "n_samples", "value": int(sample_values.nunique()) if sample_col in obs.columns else 0},
            {"metric": "sample_col", "value": sample_col if sample_col in obs.columns else ""},
            {"metric": "celltype_col", "value": celltype_col if celltype_col in obs.columns else ""},
            {"metric": "n_annotation_columns", "value": len(annotation_cols)},
            {"metric": "n_marker_positive_columns", "value": len(marker_cols)},
            {"metric": "n_ignore_annotation_columns", "value": len(ignore_cols)},
        ]
    )
    tables["overview"] = overview
    paths["overview"] = _write_table(overview, out_dir / "qc_overview.tsv")

    if sample_col in obs.columns:
        cells_per_sample = (
            sample_values.value_counts(dropna=False)
            .rename_axis(sample_col)
            .reset_index(name="n_cells")
            .sort_values(sample_col)
        )
        tables["cells_per_sample"] = cells_per_sample
        paths["cells_per_sample"] = _write_table(cells_per_sample, out_dir / "cells_per_sample.tsv")

    if celltype_col in obs.columns:
        celltype_counts = (
            obs[celltype_col].astype(str)
            .value_counts(dropna=False)
            .rename_axis(celltype_col)
            .reset_index(name="n_cells")
        )
        celltype_counts["fraction"] = celltype_counts["n_cells"] / max(int(adata.n_obs), 1)
        tables["celltype_counts"] = celltype_counts
        paths["celltype_counts"] = _write_table(celltype_counts, out_dir / "celltype_counts.tsv")

        if sample_col in obs.columns:
            by_sample = (
                obs.assign(
                    **{
                        sample_col: obs[sample_col].astype(str),
                        celltype_col: obs[celltype_col].astype(str),
                    }
                )
                .groupby([sample_col, celltype_col], observed=True)
                .size()
                .rename("n_cells")
                .reset_index()
            )
            sample_totals = sample_values.value_counts()
            by_sample["fraction_of_sample"] = by_sample["n_cells"] / by_sample[sample_col].map(sample_totals)
            tables["celltypes_per_sample"] = by_sample
            paths["celltypes_per_sample"] = _write_table(by_sample, out_dir / "celltypes_per_sample.tsv")

    if marker_cols:
        marker_rows = []
        for column in marker_cols:
            values = obs[column].fillna(False).astype(bool)
            marker_rows.append(
                {
                    "marker": str(column)[: -len(marker_suffix)],
                    "column": column,
                    "n_positive": int(values.sum()),
                    "fraction_positive": float(values.mean()),
                }
            )
        marker_table = pd.DataFrame(marker_rows)
        tables["marker_positivity"] = marker_table
        paths["marker_positivity"] = _write_table(marker_table, out_dir / "marker_positivity.tsv")

        if sample_col in obs.columns:
            rows = []
            for sample, indices in sample_values.groupby(sample_values).groups.items():
                sample_obs = obs.loc[indices]
                for column in marker_cols:
                    values = sample_obs[column].fillna(False).astype(bool)
                    rows.append(
                        {
                            sample_col: sample,
                            "marker": str(column)[: -len(marker_suffix)],
                            "column": column,
                            "n_positive": int(values.sum()),
                            "fraction_positive": float(values.mean()) if len(values) else 0.0,
                        }
                    )
            marker_by_sample = pd.DataFrame(rows)
            tables["marker_positivity_by_sample"] = marker_by_sample
            paths["marker_positivity_by_sample"] = _write_table(
                marker_by_sample,
                out_dir / "marker_positivity_by_sample.tsv",
            )

    if annotation_cols:
        rows = []
        for column in annotation_cols:
            values = obs[column].fillna(False).astype(bool)
            rows.append(
                {
                    "annotation_column": column,
                    "n_cells": int(values.sum()),
                    "fraction": float(values.mean()),
                }
            )
        annotation_table = pd.DataFrame(rows)
        tables["annotation_counts"] = annotation_table
        paths["annotation_counts"] = _write_table(annotation_table, out_dir / "annotation_counts.tsv")

    if ignore_cols:
        flagged = obs[ignore_cols].fillna(False).astype(bool).any(axis=1)
        ignore_table = pd.DataFrame(
            [
                {
                    "ignore_columns": "; ".join(map(str, ignore_cols)),
                    "n_ignore_cells": int(flagged.sum()),
                    "fraction_ignore": float(flagged.mean()),
                }
            ]
        )
        tables["ignore_cells"] = ignore_table
        paths["ignore_cells"] = _write_table(ignore_table, out_dir / "ignore_cells.tsv")

    html_parts = [
        "<html><head><meta charset=\"utf-8\"><title>QXYCell QC</title>",
        "<style>body{font-family:Arial,sans-serif;margin:24px;}table{border-collapse:collapse;margin-bottom:24px;}th,td{border:1px solid #ddd;padding:4px 8px;}th{background:#f2f2f2;}</style>",
        "</head><body><h1>QXYCell QC Report</h1>",
    ]
    for name, table in tables.items():
        html_parts.append(f"<h2>{escape(name.replace('_', ' ').title())}</h2>")
        html_parts.append(table.to_html(index=False, escape=True))
    html_parts.append("</body></html>")
    html_path = out_dir / "qc_report.html"
    html_path.write_text("\n".join(html_parts), encoding="utf-8")
    paths["html"] = html_path

    result = {
        "output_dir": out_dir,
        "paths": paths,
        "tables": tables,
    }
    adata.uns["qxycell_qc"] = {
        "output_dir": str(out_dir),
        "paths": {key: str(value) for key, value in paths.items()},
    }

    if verbose:
        print(f"Saved QXYCell QC report:\n{html_path}")

    return result
