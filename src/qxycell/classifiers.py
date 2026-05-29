"""QuPath object classifier discovery and simple-threshold parsing."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from qxycell.types import ClassifierDefinition, Message, MeasurementFile


def discover_classifier_files(project_dir: str | Path) -> list[Path]:
    """Find QuPath classifier JSON files anywhere within the project directory."""

    root = Path(project_dir).expanduser().resolve()
    files = [
        path
        for path in root.rglob("*.json")
        if path.name != "classes.json" and not path.name.startswith(".")
    ]
    return sorted(dict.fromkeys(files))


def _as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


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
                message="No QuPath object classifier JSON files were found.",
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
                message="No simple single-measurement threshold classifiers were found.",
            )
        )

    for classifier in unsupported:
        messages.append(
            Message(
                level="warning",
                code="classifiers.unsupported",
                message=(
                    f"Classifier '{classifier.path.name}' is not a simple single-measurement "
                    f"threshold classifier and will be skipped in v1 ({classifier.reason})."
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

