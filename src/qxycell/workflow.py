"""Convenience end-to-end notebook workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def workflow(
    project_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    threshold_file: str | Path | None = None,
    apply_thresholds: bool = True,
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
    """Run the common QXYCell notebook workflow and return AnnData."""

    from qxycell.celltyping import apply_celltypes
    from qxycell.filtering import remove_ignore
    from qxycell.metadata import add_metadata
    from qxycell.pipeline import apply_thresholds as apply_thresholds_to_adata, run
    from qxycell.plotting import plot_spatial, plot_stacked_bar
    from qxycell.qc import qc
    from qxycell.stage_state import checkpoint_outputs

    adata = run(
        project_dir,
        output_dir=output_dir,
        threshold_file=threshold_file,
        verbose=verbose,
    )
    if apply_thresholds:
        apply_thresholds_to_adata(
            adata,
            project_dir=project_dir,
            threshold_file=threshold_file,
            output_dir=output_dir,
            verbose=verbose,
        )

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

    adata.uns["qxycell_workflow"] = {
        "sample_metadata": str(sample_metadata) if sample_metadata is not None else None,
        "celltype_logic": str(celltype_logic) if celltype_logic is not None else None,
        "threshold_file": str(threshold_file) if threshold_file is not None else None,
        "thresholds_applied": bool(apply_thresholds),
        "remove_ignore_cells": bool(remove_ignore_cells),
        "make_qc": bool(make_qc),
        "make_plots": bool(make_plots),
    }
    checkpoint_outputs(adata)
    if verbose:
        print("Saved workflow H5AD and synchronized observation/marker tables")
    return adata
