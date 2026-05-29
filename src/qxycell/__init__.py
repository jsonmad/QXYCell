"""Python API for QXYCell."""

from qxycell.celltyping import apply_celltypes as _apply_celltypes
from qxycell.celltyping import find_latest_celltype_yaml
from qxycell.celltyping import load_celltype_logic
from qxycell.checks import CheckReport, check
from qxycell.filtering import assign_samples
from qxycell.filtering import remove_ignore
from qxycell.io_utils import load, load_latest, save
from qxycell.metadata import add_metadata
from qxycell.pipeline import run
from qxycell.plotting import plot_cn_heatmap, plot_marker_heatmap, plot_spatial, plot_stacked_bar
from qxycell.prompts import celltype_prompt
from qxycell.qc import qc
from qxycell.neighbourhoods import cn_knn, cn_kmeans, cn_name
from qxycell.tma import assign_tma_cores
from qxycell.workflow import workflow

celltype = _apply_celltypes

del _apply_celltypes

__all__ = [
    "CheckReport",
    "add_metadata",
    "assign_samples",
    "assign_tma_cores",
    "celltype",
    "cn_knn",
    "cn_kmeans",
    "cn_name",
    "celltype_prompt",
    "check",
    "find_latest_celltype_yaml",
    "load",
    "load_latest",
    "load_celltype_logic",
    "plot_cn_heatmap",
    "plot_marker_heatmap",
    "plot_spatial",
    "plot_stacked_bar",
    "qc",
    "remove_ignore",
    "run",
    "save",
    "workflow",
]
