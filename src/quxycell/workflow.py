"""Convenience end-to-end notebook workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def workflow(
    project_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    sample_metadata: str | Path | Any | None = None,
    sample_col: str = "Image",
    metadata_sample_col: str | None = None,
    celltype_logic: str | Path | dict[str, Any] | None = None,
    remove_ignore_cells: bool = True,
    make_qc: bool = True,
    make_plots: bool = True,
    plot_sample_col: str | None = None,
    verbose: bool = True,
):
    """Run the common QuXYCell notebook workflow and return AnnData."""

    from quxycell.celltyping import apply_celltypes
    from quxycell.filtering import remove_ignore
    from quxycell.io_utils import save
    from quxycell.metadata import add_metadata
    from quxycell.pipeline import run
    from quxycell.plotting import plot_spatial, plot_stacked_bar
    from quxycell.qc import qc

    adata = run(project_dir, output_dir=output_dir, verbose=verbose)

    if sample_metadata is not None:
        add_metadata(
            adata,
            sample_metadata,
            sample_col=sample_col,
            metadata_sample_col=metadata_sample_col,
            verbose=verbose,
        )

    if remove_ignore_cells:
        remove_ignore(adata, verbose=verbose)

    if celltype_logic is not None:
        apply_celltypes(adata, celltype_logic, verbose=verbose)

    if make_qc:
        qc(adata, sample_col=sample_col, verbose=verbose)

    if make_plots and "celltype" in adata.obs.columns:
        plot_col = plot_sample_col or sample_col
        plot_stacked_bar(adata, sample_col=plot_col, show=False, verbose=verbose)
        plot_spatial(adata, sample_col=plot_col, show=False, verbose=verbose)

    save(adata, verbose=verbose)
    adata.uns["quxycell_workflow"] = {
        "sample_metadata": str(sample_metadata) if sample_metadata is not None else None,
        "celltype_logic": str(celltype_logic) if celltype_logic is not None else None,
        "remove_ignore_cells": bool(remove_ignore_cells),
        "make_qc": bool(make_qc),
        "make_plots": bool(make_plots),
    }
    return adata
