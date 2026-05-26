"""Filtering helpers for QuXYCell AnnData objects."""

from __future__ import annotations


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
            print(f"No ignore annotation columns found matching '{annotation_prefix}*{ignore_text}*'.")
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
