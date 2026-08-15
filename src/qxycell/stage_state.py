"""State, invalidation, and checkpoint helpers for rerunnable QXYCell stages."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any


STAGE_ORDER = (
    "measurements",
    "annotations",
    "thresholds",
    "celltype_prompt",
    "celltypes",
    "post_analysis",
)


def _metadata(adata) -> dict[str, Any]:
    return adata.uns.setdefault("qxycell", {})


def _stages(adata) -> dict[str, Any]:
    return _metadata(adata).setdefault("stages", {})


def _remove_owned_columns(adata, record: object) -> list[str]:
    if not isinstance(record, dict):
        return []
    columns = [str(column) for column in record.get("columns", [])]
    present = [column for column in columns if column in adata.obs.columns]
    if present:
        adata.obs.drop(columns=present, inplace=True)
    return present


def _remove_owned_files(record: object) -> list[str]:
    if not isinstance(record, dict):
        return []
    removed = []
    for value in record.get("files", []):
        path = Path(str(value)).expanduser()
        if path.is_file():
            path.unlink()
            removed.append(str(path))
    return removed


def _clear_active_stage_metadata(adata, stage: str, *, remove_outputs: bool) -> None:
    metadata = _metadata(adata)
    if stage == "annotations":
        metadata.pop("annotations", None)
    elif stage == "thresholds":
        metadata["thresholding_applied"] = False
        metadata["threshold_source"] = None
        metadata["threshold_source_kind"] = "none"
        metadata["generated_threshold_template"] = None
        if remove_outputs:
            adata.uns.pop("qxycell_thresholding", None)
    elif stage == "celltype_prompt":
        metadata["llm_prompt_generated"] = False
        metadata["llm_prompt_path"] = None
    elif stage == "celltypes":
        metadata["celltyping_applied"] = False
        if remove_outputs:
            adata.uns.pop("qxycell_celltyping", None)


def prepare_stage(
    adata,
    stage: str,
    *,
    remove_downstream_columns: bool = True,
) -> dict[str, list[str]]:
    """Remove outputs owned by ``stage`` and mark all downstream stages stale."""

    if stage not in STAGE_ORDER:
        raise ValueError(f"Unknown QXYCell stage: {stage}")

    stages = _stages(adata)
    removed_columns = _remove_owned_columns(adata, stages.get(stage))
    removed_files = _remove_owned_files(stages.get(stage))
    _clear_active_stage_metadata(adata, stage, remove_outputs=True)
    stages[stage] = {
        "status": "running",
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "columns": [],
        "files": [],
    }

    stage_index = STAGE_ORDER.index(stage)
    for downstream in STAGE_ORDER[stage_index + 1 :]:
        record = stages.get(downstream)
        if not isinstance(record, dict):
            continue
        if remove_downstream_columns:
            removed_columns.extend(_remove_owned_columns(adata, record))
            removed_files.extend(_remove_owned_files(record))
        _clear_active_stage_metadata(
            adata,
            downstream,
            remove_outputs=remove_downstream_columns,
        )
        record["status"] = "stale"
        record["stale_because"] = stage
        record["stale_at"] = datetime.now().isoformat(timespec="seconds")

    return {
        "removed_columns": sorted(dict.fromkeys(removed_columns)),
        "removed_files": sorted(dict.fromkeys(removed_files)),
    }


def checkpoint_outputs(adata) -> Path:
    """Write the current AnnData and its public observation/marker tables."""

    from qxycell.io_utils import _default_h5ad_path, save
    from qxycell.paths import resolve_output_dir

    metadata = _metadata(adata)
    output_path = resolve_output_dir(adata=adata)
    tables_dir = output_path / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    metadata["output_dir"] = str(output_path)
    metadata["run_dir"] = str(output_path)
    metadata["tables_dir"] = str(tables_dir)
    h5ad_path = _default_h5ad_path(adata, output_path)
    metadata["h5ad_path"] = str(h5ad_path)
    save(adata, path=h5ad_path, verbose=False)
    adata.obs.to_csv(tables_dir / "cells_obs.csv")
    adata.var.to_csv(tables_dir / "markers_var.csv")
    return h5ad_path


def complete_stage(
    adata,
    stage: str,
    *,
    columns: list[str] | tuple[str, ...] = (),
    files: list[str | Path] | tuple[str | Path, ...] = (),
    details: dict[str, Any] | None = None,
    checkpoint: bool = True,
) -> dict[str, Any]:
    """Record successful stage ownership and optionally checkpoint the dataset."""

    if stage not in STAGE_ORDER:
        raise ValueError(f"Unknown QXYCell stage: {stage}")
    record: dict[str, Any] = {
        "status": "complete",
        "completed_at": datetime.now().isoformat(timespec="seconds"),
        "columns": sorted(dict.fromkeys(str(column) for column in columns)),
        "files": sorted(dict.fromkeys(str(Path(path).expanduser().resolve()) for path in files)),
    }
    if details:
        record.update(details)
    _stages(adata)[stage] = record
    if checkpoint:
        checkpoint_outputs(adata)
    return record
