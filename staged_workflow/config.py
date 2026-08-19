"""Edit these settings once before running the staged workflow examples."""

from pathlib import Path


PROJECT_DIR = Path("/path/to/qupath_project")
OUTPUT_DIR = Path("/path/to/qxycell_output")
THRESHOLD_TABLE = PROJECT_DIR / "thresholds.tsv"
CELLTYPE_YAML = OUTPUT_DIR / "celltype" / "celltype_logic.yaml"
PIXEL_SIZE_UM = 0.28
IGNORE_ANNOTATION_TEXT = "ignore"
CELLTYPE_CONTEXT = (
    "Describe the tissue, disease, experimental groups, and expected cell "
    "populations here. The generated rules require expert review."
)
