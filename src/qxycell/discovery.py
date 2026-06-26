"""Shared project-file discovery helpers."""

from __future__ import annotations

from pathlib import Path
import re


QXY_OUTPUT_DIR_NAMES = {"output", "outputs"}
QXY_PROJECT_OUTPUT_RE = re.compile(r".*_(?:run|check)_\d{6}_\d{4}$")


def is_qxy_output_artifact(path: str | Path, project_dir: str | Path) -> bool:
    """Return True when ``path`` is inside a generated QXYCell output folder."""

    resolved_path = Path(path).expanduser().resolve()
    resolved_root = Path(project_dir).expanduser().resolve()
    try:
        relative = resolved_path.relative_to(resolved_root)
    except ValueError:
        return False

    for part in relative.parts[:-1]:
        lower = part.lower()
        if (
            lower in QXY_OUTPUT_DIR_NAMES
            or lower.startswith("qxy_outputs_")
            or QXY_PROJECT_OUTPUT_RE.fullmatch(part)
        ):
            return True
    return False
