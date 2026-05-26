"""Python API for QuXYCell."""

from quxycell.celltyping import apply_celltypes as _apply_celltypes
from quxycell.celltyping import find_latest_celltype_yaml
from quxycell.celltyping import load_celltype_logic
from quxycell.checks import CheckReport, check
from quxycell.filtering import assign_samples
from quxycell.filtering import remove_ignore
from quxycell.io_utils import load, load_latest, save
from quxycell.metadata import add_metadata
from quxycell.pipeline import run
from quxycell.plotting import plot_spatial_celltypes, plot_stacked_bar
from quxycell.prompts import celltype_prompt
from quxycell.qc import qc
from quxycell.tma import assign_tma_cores
from quxycell.workflow import workflow

celltype = _apply_celltypes

del _apply_celltypes

__all__ = [
    "CheckReport",
    "add_metadata",
    "assign_samples",
    "assign_tma_cores",
    "celltype",
    "celltype_prompt",
    "check",
    "find_latest_celltype_yaml",
    "load",
    "load_latest",
    "load_celltype_logic",
    "plot_spatial_celltypes",
    "plot_stacked_bar",
    "qc",
    "remove_ignore",
    "run",
    "save",
    "workflow",
]
