"""QuPath object classifier discovery and simple-threshold parsing."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from qxycell.types import ClassifierDefinition, Message, MeasurementFile


MANUAL_THRESHOLD_FILENAMES = {
    "classifier_thresholds.csv",
    "classifier_thresholds.tsv",
    "manual_thresholds.csv",
    "manual_thresholds.tsv",
    "marker_thresholds.csv",
    "marker_thresholds.tsv",
    "qxycell_thresholds.csv",
    "qxycell_thresholds.tsv",
    "thresholds.csv",
    "thresholds.tsv",
}
MANUAL_THRESHOLD_PREFIX = "thresholds_"
NAME_COLUMNS = ("name", "marker", "classifier", "classifier_name")
IMAGE_COLUMNS = ("image", "Image", "sample", "Sample")
MEASUREMENT_COLUMNS = (
    "measurement_column",
    "measurement",
    "source_measurement_column",
    "measurement_col",
)
THRESHOLD_COLUMNS = ("threshold", "cutoff", "cut_off")
BASE_THRESHOLD_TEMPLATE_COLUMNS = {
    "compartment",
    "localization",
    "localisation",
    "name",
    "marker",
    "classifier",
    "classifier_name",
    "image",
    "sample",
    "measurement_column",
    "measurement",
    "source_measurement_column",
    "measurement_col",
    "threshold",
    "cutoff",
    "cut_off",
}
MEASUREMENT_COMPARTMENTS = {"Cell", "Cytoplasm", "Membrane", "Nucleus"}
MISSING_THRESHOLD_VALUES = {"", "#n/a", "n/a", "na", "nan", "none", "null", "<na>"}


def discover_classifier_files(project_dir: str | Path) -> list[Path]:
    """Find QuPath classifier JSON files anywhere within the project directory."""

    root = Path(project_dir).expanduser().resolve()
    files = [
        path
        for path in root.rglob("*.json")
        if path.name != "classes.json" and not path.name.startswith(".")
    ]
    return sorted(dict.fromkeys(files))


def discover_threshold_files(project_dir: str | Path) -> list[Path]:
    """Find manual classifier-threshold CSV/TSV files below a project directory."""

    root = Path(project_dir).expanduser().resolve()
    files = [
        path
        for path in root.rglob("*")
        if path.is_file()
        and not path.name.startswith(".")
        and _is_manual_threshold_file(path)
    ]
    return sorted(dict.fromkeys(files))


def select_threshold_file(paths: list[Path]) -> tuple[Path | None, list[Path]]:
    """Choose one manual threshold file deterministically.

    Timestamped ``thresholds_*.tsv/csv`` files are preferred. Within the chosen
    class, the most recently modified file wins. Remaining candidates are
    returned so callers can report that they were ignored.
    """

    if not paths:
        return None, []
    timestamped = [path for path in paths if _is_timestamped_threshold_file(path)]
    candidates = timestamped if timestamped else paths
    chosen = max(candidates, key=lambda path: (path.stat().st_mtime, str(path)))
    ignored = [path for path in paths if path != chosen]
    return chosen, sorted(ignored)


def _is_timestamped_threshold_file(path: Path) -> bool:
    return (
        path.suffix.lower() in {".tsv", ".csv"}
        and path.stem.lower().startswith(MANUAL_THRESHOLD_PREFIX)
    )


def _is_manual_threshold_file(path: Path) -> bool:
    name = path.name.lower()
    suffix = path.suffix.lower()
    return (
        name in MANUAL_THRESHOLD_FILENAMES
        or (
            suffix in {".tsv", ".csv"}
            and path.stem.lower().startswith(MANUAL_THRESHOLD_PREFIX)
        )
    )


def marker_name_from_measurement_column(column: str) -> str:
    """Return a readable marker name from a QuPath measurement column."""

    _compartment, text = measurement_parts_from_measurement_column(column)
    for separator in (" - ", ":"):
        if separator in text:
            text = text.split(separator, 1)[0].strip()
            break
    return text or "marker"


def measurement_parts_from_measurement_column(column: str) -> tuple[str, str]:
    """Return ``(compartment, marker_text)`` from a QuPath measurement column."""

    text = str(column).strip()
    if ":" not in text:
        return "", text

    prefix, remainder = (part.strip() for part in text.split(":", 1))
    if prefix in MEASUREMENT_COMPARTMENTS:
        return prefix, remainder
    return "", text


def compartment_from_measurement_column(column: str) -> str:
    """Return the localization compartment from a QuPath measurement column."""

    compartment, _marker = measurement_parts_from_measurement_column(column)
    return compartment


def measurement_columns_for_threshold_template(
    measurement_files: list[MeasurementFile],
) -> list[str]:
    """Return mean/median measurement columns for a manual threshold template."""

    excluded = {
        "Image",
        "Object ID",
        "Centroid X µm",
        "Centroid Y µm",
    }
    columns: list[str] = []
    for measurement_file in measurement_files:
        for column in measurement_file.columns:
            if column in excluded:
                continue
            if str(column).startswith("quxy_"):
                continue
            if not any(term in str(column).lower() for term in ("mean", "median")):
                continue
            columns.append(str(column))
    return sorted(dict.fromkeys(columns))


def _as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _is_missing_threshold_value(value: Any) -> bool:
    return str(value).strip().lower() in MISSING_THRESHOLD_VALUES


def parse_classifier(path: str | Path) -> ClassifierDefinition:
    """Parse simple QuPath ClassifyByMeasurementFunction classifiers."""

    path = Path(path).expanduser().resolve()
    try:
        data = json.loads(path.read_text(errors="replace"))
    except Exception as exc:
        return ClassifierDefinition(
            path=path,
            name=path.stem,
            measurement_column=None,
            threshold=None,
            is_simple=False,
            reason=f"unreadable_json: {exc}",
        )

    function = data.get("function") if isinstance(data, dict) else None
    if not isinstance(function, dict):
        return ClassifierDefinition(
            path=path,
            name=path.stem,
            measurement_column=None,
            threshold=None,
            is_simple=False,
            reason="missing_function",
        )

    classifier_fun = function.get("classifier_fun")
    measurement = function.get("measurement")
    threshold = _as_float(function.get("threshold"))
    if classifier_fun != "ClassifyByMeasurementFunction" or not measurement or threshold is None:
        return ClassifierDefinition(
            path=path,
            name=path.stem,
            measurement_column=measurement if isinstance(measurement, str) else None,
            threshold=threshold,
            is_simple=False,
            reason="not_single_measurement_threshold",
        )

    return ClassifierDefinition(
        path=path,
        name=path.stem,
        measurement_column=str(measurement),
        threshold=threshold,
        is_simple=True,
    )


def parse_classifiers(paths: list[Path]) -> list[ClassifierDefinition]:
    """Parse all discovered classifier files."""

    return [parse_classifier(path) for path in paths]


def _delimiter_for_threshold_path(path: Path) -> str:
    return "\t" if path.suffix.lower() == ".tsv" else ","


def _get_first(row: dict[str, str], aliases: tuple[str, ...]) -> str | None:
    lower_lookup = {str(key).strip().lower(): value for key, value in row.items()}
    for alias in aliases:
        value = lower_lookup.get(alias.lower())
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def _row_has_threshold_value(row: dict[str, str]) -> bool:
    threshold_value = _get_first(row, THRESHOLD_COLUMNS)
    if threshold_value and not _is_missing_threshold_value(threshold_value):
        return True
    for column, value in row.items():
        column_name = str(column).strip()
        if column_name.lower() in BASE_THRESHOLD_TEMPLATE_COLUMNS:
            continue
        if str(value).strip() and not _is_missing_threshold_value(value):
            return True
    return False


def parse_threshold_file(path: str | Path) -> list[ClassifierDefinition]:
    """Parse manually entered marker threshold CSV/TSV rows.

    Required columns are a marker/name column, a measurement column, and a
    numeric threshold column. Accepted aliases are defined by ``NAME_COLUMNS``,
    ``MEASUREMENT_COLUMNS``, and ``THRESHOLD_COLUMNS``.
    """

    path = Path(path).expanduser().resolve()
    delimiter = _delimiter_for_threshold_path(path)
    classifiers: list[ClassifierDefinition] = []

    try:
        with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row_index, row in enumerate(reader, start=2):
                name = _get_first(row, NAME_COLUMNS)
                measurement = _get_first(row, MEASUREMENT_COLUMNS)
                row_has_threshold = _row_has_threshold_value(row)
                if not row_has_threshold:
                    continue
                if not name or not measurement:
                    classifiers.append(
                        ClassifierDefinition(
                            path=path.parent / f"{path.name}#row{row_index}",
                            name=name or f"{path.stem}_row_{row_index}",
                            measurement_column=measurement,
                            threshold=None,
                            is_simple=False,
                            reason="invalid_manual_threshold_row",
                        )
                    )
                    continue

                image_thresholds = _wide_image_thresholds(row)
                if image_thresholds:
                    for image, threshold in image_thresholds:
                        classifiers.append(
                            ClassifierDefinition(
                                path=path.parent / f"{path.name}#row{row_index}:{image}",
                                name=name,
                                measurement_column=measurement,
                                threshold=threshold,
                                is_simple=True,
                                image=image,
                            )
                        )
                    continue

                image = _get_first(row, IMAGE_COLUMNS)
                threshold = _as_float(_get_first(row, THRESHOLD_COLUMNS))
                row_path = path.parent / f"{path.name}#row{row_index}"
                if threshold is None:
                    classifiers.append(
                        ClassifierDefinition(
                            path=row_path,
                            name=name,
                            measurement_column=measurement,
                            threshold=threshold,
                            is_simple=False,
                            image=image,
                            reason="invalid_manual_threshold_row",
                        )
                    )
                    continue
                classifiers.append(
                    ClassifierDefinition(
                        path=row_path,
                        name=name,
                        measurement_column=measurement,
                        threshold=threshold,
                        is_simple=True,
                        image=image,
                    )
                )
    except Exception as exc:
        return [
            ClassifierDefinition(
                path=path,
                name=path.stem,
                measurement_column=None,
                threshold=None,
                is_simple=False,
                reason=f"unreadable_manual_threshold_file: {exc}",
            )
        ]

    return classifiers


def _wide_image_thresholds(row: dict[str, str]) -> list[tuple[str, float]]:
    """Return image-specific threshold values from wide manual-template columns."""

    thresholds: list[tuple[str, float]] = []
    for column, value in row.items():
        column_name = str(column).strip()
        if column_name.lower() in BASE_THRESHOLD_TEMPLATE_COLUMNS:
            continue
        threshold = _as_float(value)
        if threshold is None:
            continue
        thresholds.append((column_name, threshold))
    return thresholds


def parse_threshold_files(paths: list[Path]) -> list[ClassifierDefinition]:
    """Parse all manual threshold files."""

    classifiers: list[ClassifierDefinition] = []
    for path in paths:
        classifiers.extend(parse_threshold_file(path))
    return classifiers


def validate_classifiers(
    classifiers: list[ClassifierDefinition],
    measurement_files: list[MeasurementFile],
) -> list[Message]:
    """Validate classifier references against measurement columns."""

    messages: list[Message] = []
    if not classifiers:
        messages.append(
            Message(
                level="error",
                code="classifiers.missing",
                message="No classifier JSON files or manual threshold TSV/CSV files were found.",
            )
        )
        return messages

    simple = [classifier for classifier in classifiers if classifier.is_simple]
    unsupported = [classifier for classifier in classifiers if not classifier.is_simple]
    if not simple:
        messages.append(
            Message(
                level="error",
                code="classifiers.no_simple_thresholds",
                message="No usable single-measurement threshold classifiers were found.",
            )
        )

    for classifier in unsupported:
        messages.append(
            Message(
                level="warning",
                code="classifiers.unsupported",
                message=(
                    f"Classifier '{classifier.path.name}' is not a simple single-measurement "
                    f"threshold definition and will be skipped in v1 ({classifier.reason})."
                ),
                path=str(classifier.path),
            )
        )

    all_columns = set()
    for measurement_file in measurement_files:
        all_columns.update(measurement_file.columns)

    for classifier in simple:
        if classifier.measurement_column not in all_columns:
            messages.append(
                Message(
                    level="error",
                    code="classifiers.measurement_column_missing",
                    message=(
                        f"Classifier '{classifier.path.name}' references missing measurement "
                        f"column: {classifier.measurement_column}"
                    ),
                    path=str(classifier.path),
                )
            )
    return messages
