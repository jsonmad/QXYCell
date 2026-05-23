"""Python API for QUXYCell."""

from quxycell.celltyping import apply_celltypes as _apply_celltypes
from quxycell.celltyping import find_latest_celltype_yaml
from quxycell.celltyping import load_celltype_logic
from quxycell.checks import CheckReport, check
from quxycell.filtering import remove_ignore
from quxycell.pipeline import run
from quxycell.plotting import plot_spatial_celltypes, plot_stacked_bar
from quxycell.prompts import celltype_prompt

celltype = _apply_celltypes

del _apply_celltypes

__all__ = [
    "CheckReport",
    "celltype",
    "celltype_prompt",
    "check",
    "find_latest_celltype_yaml",
    "load_celltype_logic",
    "plot_spatial_celltypes",
    "plot_stacked_bar",
    "remove_ignore",
    "run",
]
