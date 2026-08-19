"""Python API for QXYCell."""

from qxycell.celltyping import apply_celltypes
from qxycell.celltyping import find_latest_celltype_yaml
from qxycell.celltyping import load_celltype_logic
from qxycell.checks import CheckReport, check, generate_threshold_table
from qxycell.geojson import load_cell_polygons
from qxycell.filtering import assign_annotations
from qxycell.filtering import assign_core_ids_from_measurements
from qxycell.filtering import assign_samples
from qxycell.filtering import remove_annotations
from qxycell.filtering import remove_cells
from qxycell.io_utils import load, load_latest, save
from qxycell.metadata import add_metadata
from qxycell.pipeline import (
    add_annotations,
    apply_thresholds,
    import_cells,
    threshold,
    threshold_from_classifiers,
    threshold_from_table,
)
from qxycell.plotting import (
    plot_annotation_polygons,
    plot_cell_boundaries,
    plot_cn_heatmap,
    plot_marker_intensity_heatmap,
    plot_marker_positivity_heatmap,
    plot_marker_heatmap,
    plot_spatial,
    plot_stacked_bar,
)
from qxycell.prompts import celltype_prompt
from qxycell.qc import qc
from qxycell.neighbourhoods import cn_knn, cn_kmeans, cn_name
from qxycell.workflow import workflow

celltype = apply_celltypes

__all__ = [
    "CheckReport",
    "add_metadata",
    "add_annotations",
    "apply_celltypes",
    "apply_thresholds",
    "assign_annotations",
    "assign_core_ids_from_measurements",
    "assign_samples",
    "celltype",
    "cn_knn",
    "cn_kmeans",
    "cn_name",
    "celltype_prompt",
    "check",
    "find_latest_celltype_yaml",
    "generate_threshold_table",
    "import_cells",
    "load",
    "load_cell_polygons",
    "load_latest",
    "load_celltype_logic",
    "plot_cn_heatmap",
    "plot_annotation_polygons",
    "plot_cell_boundaries",
    "plot_marker_heatmap",
    "plot_marker_intensity_heatmap",
    "plot_marker_positivity_heatmap",
    "plot_spatial",
    "plot_stacked_bar",
    "qc",
    "remove_annotations",
    "remove_cells",
    "save",
    "threshold",
    "threshold_from_classifiers",
    "threshold_from_table",
    "workflow",
]
