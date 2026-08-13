"""Shared path helpers for QXYCell outputs."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import re


OUTPUT_DIR_PREFIX = "qxy_outputs_"
OUTPUT_TIMESTAMP_FORMAT = "%y%m%d-%H%M"
PROJECT_OUTPUT_TIMESTAMP_FORMAT = "%y%m%d_%H%M"
LEGACY_DEFAULT_OUTPUT_PARTS = {
    ("outputs", "qxy" "_check"),
    ("outputs", "qxy" "_run"),
}


def _is_legacy_default_output_path(path: Path) -> bool:
    return any(
        path.parts[-len(parts):] == parts
        for parts in LEGACY_DEFAULT_OUTPUT_PARTS
    )


def timestamped_output_dir(base_dir: str | Path = ".") -> Path:
    """Return the default timestamped QXYCell output directory."""

    timestamp = datetime.now().strftime(OUTPUT_TIMESTAMP_FORMAT)
    return (Path(base_dir).expanduser() / f"{OUTPUT_DIR_PREFIX}{timestamp}").resolve()


def timestamped_project_output_dir(project_dir: str | Path, *, kind: str = "run") -> Path:
    """Return the default output directory beside a QuPath project folder."""

    project_path = Path(project_dir).expanduser().resolve()
    timestamp = datetime.now().strftime(PROJECT_OUTPUT_TIMESTAMP_FORMAT)
    safe_kind = str(kind).strip().lower() or "run"
    return (project_path.parent / f"{project_path.name}_{safe_kind}_{timestamp}").resolve()


def latest_timestamped_output_dir(base_dir: str | Path = ".") -> Path | None:
    """Return the most recently modified timestamped QXYCell run directory."""

    base_path = Path(base_dir).expanduser().resolve()
    candidates = [
        path
        for path in list(base_path.glob(f"{OUTPUT_DIR_PREFIX}*")) + list(base_path.glob("*_run_*"))
        if path.is_dir()
        and (
            path.name.startswith(OUTPUT_DIR_PREFIX)
            or re.fullmatch(r".*_run_\d{6}_\d{4}", path.name)
        )
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def output_dir_from_adata(adata) -> Path | None:
    """Return the QXYCell output directory stored on an AnnData object."""

    metadata = getattr(adata, "uns", {}).get("qxycell", {})
    output_dir = metadata.get("output_dir") if isinstance(metadata, dict) else None
    if not output_dir:
        return None
    output_path = Path(output_dir).expanduser().resolve()
    if _is_legacy_default_output_path(output_path):
        return None
    return output_path


def resolve_output_dir(
    output_dir: str | Path | None = None,
    *,
    adata=None,
    project_dir: str | Path | None = None,
    project_output_kind: str = "run",
) -> Path:
    """Resolve an explicit output path or the default timestamped output path."""

    if output_dir is not None:
        output_path = Path(output_dir).expanduser().resolve()
        if _is_legacy_default_output_path(output_path):
            return (
                timestamped_project_output_dir(
                    project_dir,
                    kind=project_output_kind,
                )
                if project_dir is not None
                else timestamped_output_dir()
            )
        return output_path
    if adata is not None:
        adata_output_dir = output_dir_from_adata(adata)
        if adata_output_dir is not None:
            return adata_output_dir
    if project_dir is not None:
        return timestamped_project_output_dir(project_dir, kind=project_output_kind)
    return timestamped_output_dir()
