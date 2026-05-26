"""Filtering helpers for QuXYCell AnnData objects."""

from __future__ import annotations

import json
import warnings


def _annotation_label_for_column(adata, column: str, annotation_prefix: str) -> str:
    label_maps = [
        getattr(adata, "uns", {}).get("quxycell_annotation_labels", {}),
        getattr(adata, "uns", {}).get("quxycell", {}).get("annotation_labels", {}),
    ]
    for label_map in label_maps:
        if isinstance(label_map, dict) and column in label_map:
            return str(label_map[column])

    label = str(column)
    if label.startswith(annotation_prefix):
        label = label[len(annotation_prefix):]
    return label


def assign_samples(
    adata,
    *,
    annotation_prefix: str = "annotation__",
    sample_text: str = "sample",
    sample_col: str = "Sample",
    unassigned_label: str | None = "Unassigned",
    conflict_label: str = "Ambiguous",
    overwrite: bool = True,
    verbose: bool = True,
    max_conflict_examples: int = 10,
) -> dict[str, object]:
    """Assign one sample label per cell from sample annotation columns.

    The function searches ``adata.obs`` for annotation columns whose names start
    with ``annotation_prefix`` and contain ``sample_text`` case-insensitively.
    Cells inside exactly one matching annotation receive that annotation label
    in ``adata.obs[sample_col]``. Cells inside multiple sample annotations are
    labelled with ``conflict_label`` and emit a warning.
    """

    if sample_col in adata.obs.columns and not overwrite:
        raise ValueError(
            f"adata.obs already contains {sample_col!r}. "
            "Pass overwrite=True or choose a different sample_col."
        )

    sample_columns = [
        column
        for column in adata.obs.columns
        if str(column).startswith(annotation_prefix)
        and sample_text.lower() in str(column).lower()
    ]
    sample_columns = sorted(sample_columns)

    if not sample_columns:
        adata.obs[sample_col] = unassigned_label
        summary = {
            "sample_col": sample_col,
            "sample_text": sample_text,
            "annotation_columns": [],
            "sample_names": [],
            "n_assigned_cells": 0,
            "n_unassigned_cells": int(adata.n_obs),
            "n_conflicting_cells": 0,
            "conflict_examples": [],
        }
        adata.uns["quxycell_sample_annotations"] = {
            **summary,
            "conflict_examples": json.dumps([]),
        }
        if verbose:
            print(
                "No sample annotation columns found matching "
                f"'{annotation_prefix}*{sample_text}*'."
            )
            print(f"Created {sample_col!r} with {unassigned_label!r}.")
        return summary

    hits = adata.obs[sample_columns].fillna(False).astype(bool)
    hit_counts = hits.sum(axis=1)
    sample_names = {
        column: _annotation_label_for_column(adata, str(column), annotation_prefix)
        for column in sample_columns
    }

    assigned = adata.obs.iloc[:, 0].astype(object).copy()
    assigned[:] = unassigned_label
    for column in sample_columns:
        assigned.loc[hits[column] & (hit_counts == 1)] = sample_names[column]

    conflict_mask = hit_counts > 1
    if bool(conflict_mask.any()):
        assigned.loc[conflict_mask] = conflict_label

    adata.obs[sample_col] = assigned

    conflict_examples = []
    if bool(conflict_mask.any()):
        for index in hits.index[conflict_mask][:max_conflict_examples]:
            columns = [column for column in sample_columns if bool(hits.at[index, column])]
            conflict_examples.append(
                {
                    "obs_index": str(index),
                    "sample_annotations": [sample_names[column] for column in columns],
                    "annotation_columns": [str(column) for column in columns],
                }
            )
        warnings.warn(
            f"{int(conflict_mask.sum()):,} cells fall inside more than one sample annotation. "
            f"They were labelled {conflict_label!r} in adata.obs[{sample_col!r}].",
            stacklevel=2,
        )

    summary = {
        "sample_col": sample_col,
        "sample_text": sample_text,
        "annotation_columns": [str(column) for column in sample_columns],
        "sample_names": [sample_names[column] for column in sample_columns],
        "n_assigned_cells": int((hit_counts == 1).sum()),
        "n_unassigned_cells": int((hit_counts == 0).sum()),
        "n_conflicting_cells": int(conflict_mask.sum()),
        "conflict_examples": conflict_examples,
    }
    adata.uns["quxycell_sample_annotations"] = {
        **summary,
        "conflict_examples": json.dumps(conflict_examples),
    }

    if verbose:
        print(f"Assigned sample labels to adata.obs[{sample_col!r}]")
        print(f"Sample annotation columns: {len(sample_columns)}")
        print(f"Assigned cells: {summary['n_assigned_cells']:,}")
        print(f"Unassigned cells: {summary['n_unassigned_cells']:,}")
        print(f"Conflicting cells: {summary['n_conflicting_cells']:,}")

    return summary


def remove_ignore(
    adata,
    *,
    annotation_prefix: str = "annotation__",
    ignore_text: str = "Ignore",
    copy: bool = False,
    verbose: bool = True,
):
    """Remove cells that fall inside annotation columns containing ``Ignore``.

    The function searches ``adata.obs`` for annotation columns whose names start
    with ``annotation_prefix`` and contain ``ignore_text`` case-insensitively.
    A cell is removed if any matching column is true. By default, the input
    AnnData object is modified in place and returned. Use ``copy=True`` to
    return a filtered copy without changing the original object.
    """

    columns = [
        column
        for column in adata.obs.columns
        if str(column).startswith(annotation_prefix)
        and ignore_text.lower() in str(column).lower()
    ]
    if not columns:
        if verbose:
            print(
                "No ignore annotation columns found matching "
                f"'{annotation_prefix}*{ignore_text}*'."
            )
        return adata.copy() if copy else adata

    flagged = adata.obs[columns].fillna(False).astype(bool).any(axis=1)
    n_removed = int(flagged.sum())
    n_before = int(adata.n_obs)
    keep = ~flagged.to_numpy()
    if copy:
        filtered = adata[keep].copy()
    else:
        adata._inplace_subset_obs(keep)
        filtered = adata

    if verbose:
        print("Removed cells in ignore annotations")
        print(f"Ignore columns: {columns}")
        print(f"Cells before: {n_before:,}")
        print(f"Cells removed: {n_removed:,}")
        print(f"Cells after: {filtered.n_obs:,}")

    return filtered
