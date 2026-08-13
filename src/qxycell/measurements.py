"""Measurement table discovery and validation."""

from __future__ import annotations

import csv
from pathlib import Path

from qxycell.discovery import is_qxy_output_artifact
from qxycell.types import MeasurementFile, Message

REQUIRED_MEASUREMENT_COLUMNS = ("Image", "Object ID", "Centroid X µm", "Centroid Y µm")
MEASUREMENT_COLUMN_ALIASES = {
    "Centroid X ¬µm": "Centroid X µm",
    "Centroid Y ¬µm": "Centroid Y µm",
}
MEASUREMENT_TEXT_ENCODING = "utf-8-sig"


def discover_measurement_files(project_dir: str | Path) -> list[Path]:
    """Find likely QuPath measurement tables below a project export directory."""

    root = Path(project_dir).expanduser().resolve()
    candidates: list[Path] = []
    for suffix in ("*.tsv", "*.csv"):
        candidates.extend(root.rglob(suffix))

    filtered = []
    for path in candidates:
        lower = path.name.lower()
        if lower.startswith("."):
            continue
        if is_qxy_output_artifact(path, root):
            continue
        if "measurement" in lower or lower in {"detections.tsv", "detections.csv"}:
            filtered.append(path)
    return sorted(dict.fromkeys(filtered))


def required_columns() -> tuple[str, ...]:
    """Return QXYCell's required v1 QuPath measurement columns."""

    return REQUIRED_MEASUREMENT_COLUMNS


def normalize_measurement_columns(columns) -> tuple[str, ...]:
    """Normalize known QuPath measurement-header encoding variants."""

    return tuple(MEASUREMENT_COLUMN_ALIASES.get(str(column), str(column)) for column in columns)


def delimiter_for_path(path: Path) -> str:
    """Return the expected delimiter for a CSV/TSV path."""

    return "\t" if path.suffix.lower() == ".tsv" else ","


def summarize_measurement_file(path: str | Path, count_rows: bool = False) -> MeasurementFile:
    """Read a measurement table header and optionally count rows."""

    path = Path(path).expanduser().resolve()
    delimiter = delimiter_for_path(path)
    with path.open(
        newline="",
        errors="replace",
        encoding=MEASUREMENT_TEXT_ENCODING,
    ) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            columns = normalize_measurement_columns(next(reader))
        except StopIteration as exc:
            raise ValueError(f"Measurement table is empty: {path}") from exc
        n_rows = sum(1 for _ in reader) if count_rows else None
    return MeasurementFile(
        path=path,
        delimiter=delimiter,
        n_columns=len(columns),
        columns=columns,
        n_rows=n_rows,
    )


def validate_measurement_files(files: list[MeasurementFile]) -> list[Message]:
    """Validate required QuPath measurement columns."""

    messages: list[Message] = []
    if not files:
        return [
            Message(
                level="error",
                code="measurements.missing",
                message="No measurement CSV/TSV files were found.",
            )
        ]

    for file in files:
        missing = [col for col in REQUIRED_MEASUREMENT_COLUMNS if col not in file.columns]
        if missing:
            messages.append(
                Message(
                    level="error",
                    code="measurements.required_columns_missing",
                    message=(
                        "Measurement file is missing required QuPath column(s): "
                        + ", ".join(missing)
                    ),
                    path=str(file.path),
                )
            )
    return messages
