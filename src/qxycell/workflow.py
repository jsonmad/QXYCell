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
    make_dataset_summary: bool = True,
    make_plots: bool = True,
    plot_sample_col: str | None = None,
    verbose: bool = True,
):
    """Run the common QXYCell notebook workflow and return AnnData."""

    from qxycell.celltyping import apply_celltypes
    from qxycell.filtering import remove_cells
    from qxycell.metadata import add_metadata
    from qxycell.pipeline import add_annotations, import_cells
    from qxycell.pipeline import apply_thresholds as apply_thresholds_to_adata
    from qxycell.paths import resolve_output_dir
    from qxycell.plotting import plot_spatial, plot_stacked_bar
    from qxycell.summary import dataset_summary
    from qxycell.stage_state import checkpoint_outputs

    output_path = resolve_output_dir(
        output_dir,
        project_dir=project_dir,
        project_output_kind="run",
    )
    adata = import_cells(
        project_dir,
        output_dir=output_path,
        verbose=verbose,
    )
    add_annotations(adata, project_dir=project_dir, pixel_size_um=0.28, verbose=verbose)
    if apply_thresholds:
        apply_thresholds_to_adata(
            adata,
            project_dir=project_dir,
            threshold_file=threshold_file,
            output_dir=output_path,
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
        remove_cells(adata, verbose=verbose)

    if celltype_logic is not None:
        apply_celltypes(adata, celltype_logic, verbose=verbose)

    if make_dataset_summary:
        dataset_summary(
            adata,
            sample_col=sample_col,
            output_dir=output_path,
            verbose=verbose,
        )

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
        "make_dataset_summary": bool(make_dataset_summary),
        "make_plots": bool(make_plots),
    }
    checkpoint_outputs(adata)
    if verbose:
        print("Saved workflow H5AD and synchronized observation/marker tables")
    return adata
