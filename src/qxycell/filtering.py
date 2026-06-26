"""Filtering helpers for QXYCell AnnData objects."""

from __future__ import annotations

import json
import re
import warnings

import pandas as pd


def _annotation_label_for_column(adata, column: str, annotation_prefix: str) -> str:
    label_maps = [
        getattr(adata, "uns", {}).get("qxycell_annotation_labels", {}),
        getattr(adata, "uns", {}).get("qxycell", {}).get("annotation_labels", {}),
    ]
    for label_map in label_maps:
        if isinstance(label_map, dict) and column in label_map:
            return str(label_map[column])

    label = str(column)
    if label.startswith(annotation_prefix):
        label = label[len(annotation_prefix):]
    return label


def _safe_annotation_column(annotation_prefix: str, label: str) -> str:
    safe = "".join(ch if ch.isalnum() else "_" for ch in str(label).strip())
    safe = "_".join(part for part in safe.split("_") if part)
    return f"{annotation_prefix}{safe or 'Unclassified'}"


def _strip_label_count(label: str) -> str:
    return re.sub(r"\s*\(\d+\)\s*$", "", str(label).strip())


def _valid_measurement_core_values(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip()
    invalid = text.isna() | text.str.lower().isin(
        {"", "nan", "none", "<na>", "na", "null", "unassigned", "root object (image)"}
    )
    return text.mask(invalid, pd.NA)


def assign_core_ids_from_measurements(
    adata,
    *,
    target_col: str = "CoreID",
    source_cols: tuple[str, ...] = ("TMA Core", "Parent"),
    overwrite: bool = True,
    fill_only_missing: bool = False,
    drop_source_cols: bool = False,
    verbose: bool = True,
) -> dict[str, object]:
    """Assign TMA core IDs from QuPath measurement metadata columns.

    QXYCell uses this as the CoreID path when QuPath exported core metadata
    such as ``"TMA Core"`` or ``"Parent"``.
    """

    if not target_col:
        raise ValueError("target_col must be a non-empty string.")
    if isinstance(source_cols, str):
        source_cols = (source_cols,)
    source_cols = tuple(str(column) for column in source_cols if str(column).strip())
    if not source_cols:
        raise ValueError("source_cols must contain at least one column name.")

    available_source_cols = [column for column in source_cols if column in adata.obs.columns]
    missing_source_cols = [column for column in source_cols if column not in adata.obs.columns]
    if not available_source_cols:
        raise KeyError(
            "None of the requested measurement core columns exist in adata.obs: "
            f"{list(source_cols)}"
        )
    if target_col in adata.obs.columns and not overwrite and not fill_only_missing:
        raise ValueError(
            f"adata.obs already contains {target_col!r}. "
            "Pass overwrite=True, fill_only_missing=True, or choose a different target_col."
        )

    assigned = pd.Series(pd.NA, index=adata.obs.index, dtype="object")
    source_used = pd.Series(pd.NA, index=adata.obs.index, dtype="object")
    for column in available_source_cols:
        values = _valid_measurement_core_values(adata.obs[column])
        mask = assigned.isna() & values.notna()
        assigned.loc[mask] = values.loc[mask].astype(object)
        source_used.loc[mask] = column

    if fill_only_missing and target_col in adata.obs.columns:
        existing = _valid_measurement_core_values(adata.obs[target_col])
        missing_target = existing.isna()
        final_values = existing.astype(object)
        final_values.loc[missing_target] = assigned.loc[missing_target]
    else:
        final_values = assigned

    adata.obs[target_col] = pd.Categorical(final_values)

    dropped_source_cols: list[str] = []
    if drop_source_cols:
        dropped_source_cols = [
            column for column in available_source_cols if column != target_col
        ]
        if dropped_source_cols:
            adata.obs.drop(columns=dropped_source_cols, inplace=True)

    summary = {
        "target_col": target_col,
        "source_cols": list(source_cols),
        "available_source_cols": available_source_cols,
        "missing_source_cols": missing_source_cols,
        "n_assigned_cells": int(pd.Series(final_values).notna().sum()),
        "n_unassigned_cells": int(pd.Series(final_values).isna().sum()),
        "source_counts": source_used.value_counts(dropna=True).to_dict(),
        "fill_only_missing": bool(fill_only_missing),
        "dropped_source_cols": dropped_source_cols,
    }
    adata.uns["qxycell_core_ids_from_measurements"] = summary

    if verbose:
        print(f"Assigned measurement core IDs to adata.obs[{target_col!r}]")
        print(f"Source columns used: {available_source_cols}")
        if missing_source_cols:
            print(f"Source columns not present: {missing_source_cols}")
        print(f"Assigned cells: {summary['n_assigned_cells']:,}")
        print(f"Unassigned cells (NaN): {summary['n_unassigned_cells']:,}")
        if dropped_source_cols:
            print(f"Dropped source columns: {len(dropped_source_cols)}")

    return summary


def assign_annotations(
    adata,
    labels,
    *,
    target_col: str,
    annotation_prefix: str = "annotation__",
    unassigned_label: str | None = None,
    conflict_label: str = "Ambiguous",
    overwrite: bool = True,
    drop: bool = True,
    strip_counts: bool = True,
    verbose: bool = True,
    max_conflict_examples: int = 10,
) -> dict[str, object]:
    """Assign selected annotation labels to one ``adata.obs`` column.

    ``labels`` are matched against stored QXYCell annotation labels and their
    generated ``annotation__...`` column names. Labels such as ``"Tumor (26)"``
    are treated as ``"Tumor"`` by default. Cells inside exactly one selected
    annotation receive that label in ``adata.obs[target_col]``. Cells inside
    multiple selected annotations receive ``conflict_label``.
    """

    if not target_col:
        raise ValueError("target_col must be a non-empty string.")
    if target_col in adata.obs.columns and not overwrite:
        raise ValueError(
            f"adata.obs already contains {target_col!r}. "
            "Pass overwrite=True or choose a different target_col."
        )
    if isinstance(labels, str):
        labels = [label for label in labels.split(",") if label.strip()]

    requested_labels = [
        _strip_label_count(label) if strip_counts else str(label).strip()
        for label in labels
        if str(label).strip()
    ]
    if not requested_labels:
        raise ValueError("labels must contain at least one non-empty label.")

    label_map = getattr(adata, "uns", {}).get("qxycell_annotation_labels", {})
    label_to_column = {
        str(label): str(column)
        for column, label in label_map.items()
        if str(column) in adata.obs.columns
    } if isinstance(label_map, dict) else {}

    annotation_columns = {}
    missing_labels = []
    for label in requested_labels:
        candidates = [
            label_to_column.get(label),
            _safe_annotation_column(annotation_prefix, label),
            label,
        ]
        column = next(
            (candidate for candidate in candidates if candidate in adata.obs.columns),
            None,
        )
        if column is None:
            missing_labels.append(label)
        else:
            annotation_columns[label] = column

    if missing_labels:
        raise KeyError(f"Missing expected annotation columns: {missing_labels}")

    columns = list(annotation_columns.values())
    hits = adata.obs[columns].fillna(False).astype(bool)
    hit_counts = hits.sum(axis=1)

    assigned = adata.obs.iloc[:, 0].astype(object).copy()
    assigned[:] = unassigned_label
    for label, column in annotation_columns.items():
        assigned.loc[hits[column] & (hit_counts == 1)] = label

    conflict_mask = hit_counts > 1
    if bool(conflict_mask.any()):
        assigned.loc[conflict_mask] = conflict_label

    adata.obs[target_col] = pd.Categorical(assigned)

    conflict_examples = []
    if bool(conflict_mask.any()):
        for index in hits.index[conflict_mask][:max_conflict_examples]:
            hit_columns = [column for column in columns if bool(hits.at[index, column])]
            conflict_examples.append(
                {
                    "obs_index": str(index),
                    "annotation_labels": [
                        _annotation_label_for_column(adata, str(column), annotation_prefix)
                        for column in hit_columns
                    ],
                    "annotation_columns": [str(column) for column in hit_columns],
                }
            )
        warnings.warn(
            f"{int(conflict_mask.sum()):,} cells fall inside more than one selected "
            f"annotation. They were labelled {conflict_label!r} in "
            f"adata.obs[{target_col!r}].",
            stacklevel=2,
        )

    if drop:
        adata.obs.drop(columns=columns, inplace=True)
        if isinstance(label_map, dict):
            for column in columns:
                label_map.pop(column, None)
            adata.uns["qxycell_annotation_labels"] = label_map

    summary = {
        "target_col": target_col,
        "annotation_columns": [str(column) for column in columns],
        "annotation_labels": list(annotation_columns.keys()),
        "n_assigned_cells": int((hit_counts == 1).sum()),
        "n_unassigned_cells": int((hit_counts == 0).sum()),
        "n_conflicting_cells": int(conflict_mask.sum()),
        "dropped_annotation_columns": bool(drop),
        "conflict_examples": conflict_examples,
    }
    adata.uns.setdefault("qxycell_assigned_annotations", {})[target_col] = {
        **summary,
        "conflict_examples": json.dumps(conflict_examples),
    }

    if verbose:
        print(f"Assigned annotation labels to adata.obs[{target_col!r}]")
        print(f"Annotation columns: {len(columns)}")
        print(f"Assigned cells: {summary['n_assigned_cells']:,}")
        print(f"Unassigned cells (NaN): {summary['n_unassigned_cells']:,}")
        print(f"Conflicting cells ({conflict_label!r}): {summary['n_conflicting_cells']:,}")
        if drop:
            print(f"Dropped annotation columns: {len(columns)}")

    return summary


def assign_samples(
    adata,
    *,
    annotation_prefix: str = "annotation__",
    sample_text: str = "sample",
    sample_col: str = "Sample",
    unassigned_label: str | None = None,
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
        adata.obs[sample_col] = pd.Categorical([unassigned_label] * adata.n_obs)
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
        adata.uns["qxycell_sample_annotations"] = {
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

    # Convert to Categorical so unassigned cells (NaN) are proper missing data
    # and do not appear in value_counts, groupby, or cat.categories.
    adata.obs[sample_col] = pd.Categorical(assigned)

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
    adata.uns["qxycell_sample_annotations"] = {
        **summary,
        "conflict_examples": json.dumps(conflict_examples),
    }

    if verbose:
        print(f"Assigned sample labels to adata.obs[{sample_col!r}]")
        print(f"Sample annotation columns: {len(sample_columns)}")
        print(f"Assigned cells: {summary['n_assigned_cells']:,}")
        print(f"Unassigned cells (NaN): {summary['n_unassigned_cells']:,}")
        print(f"Conflicting cells ({conflict_label!r}): {summary['n_conflicting_cells']:,}")

    return summary


def remove_annotations(
    adata,
    *,
    annotation_prefix: str = "annotation__",
    text: str = "Ignore",
    copy: bool = False,
    verbose: bool = True,
):
    """Remove cells that fall inside annotation columns containing ``text``.

    The function searches ``adata.obs`` for annotation columns whose names start
    with ``annotation_prefix`` and contain ``text`` case-insensitively.
    A cell is removed if any matching column is true. By default, the input
    AnnData object is modified in place and returned. Use ``copy=True`` to
    return a filtered copy without changing the original object.
    """

    if not text:
        raise ValueError("text must be a non-empty string.")

    columns = [
        column
        for column in adata.obs.columns
        if str(column).startswith(annotation_prefix)
        and text.lower() in str(column).lower()
    ]
    if not columns:
        if verbose:
            print(
                "No annotation columns found matching "
                f"'{annotation_prefix}*{text}*'."
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
        print("Removed cells in matching annotations")
        print(f"Annotation text: {text!r}")
        print(f"Annotation columns: {columns}")
        print(f"Cells before: {n_before:,}")
        print(f"Cells removed: {n_removed:,}")
        print(f"Cells after: {filtered.n_obs:,}")

    return filtered


def remove_ignore(
    adata,
    *,
    annotation_prefix: str = "annotation__",
    ignore_text: str = "Ignore",
    copy: bool = False,
    verbose: bool = True,
):
    """Remove cells that fall inside annotation columns containing ``ignore_text``.

    This is a convenience wrapper around :func:`remove_annotations` using
    ``ignore_text="Ignore"`` by default.
    """

    return remove_annotations(
        adata,
        annotation_prefix=annotation_prefix,
        text=ignore_text,
        copy=copy,
        verbose=verbose,
    )
