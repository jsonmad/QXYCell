"""Project check/preflight workflow."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from quxycell.classifiers import (
    discover_classifier_files,
    parse_classifiers,
    validate_classifiers,
)
from quxycell.geojson import discover_geojson_files, summarize_geojson_files, validate_geojson_files
from quxycell.measurements import (
    discover_measurement_files,
    summarize_measurement_file,
    validate_measurement_files,
)
from quxycell.types import ClassifierDefinition, GeoJsonFile, MeasurementFile, Message


@dataclass(frozen=True)
class CheckReport:
    """Combined inventory and validation report for a QuPath export."""

    project_dir: Path
    output_dir: Path
    measurement_files: list[MeasurementFile]
    classifiers: list[ClassifierDefinition]
    geojson_files: list[GeoJsonFile]
    messages: list[Message]

    @property
    def ok(self) -> bool:
        """True when no errors were emitted."""

        return not any(message.level == "error" for message in self.messages)

    @property
    def n_errors(self) -> int:
        return sum(1 for message in self.messages if message.level == "error")

    @property
    def n_warnings(self) -> int:
        return sum(1 for message in self.messages if message.level == "warning")

    def to_dict(self) -> dict[str, Any]:
        """Serialize report content."""

        return {
            "project_dir": str(self.project_dir),
            "output_dir": str(self.output_dir),
            "ok": self.ok,
            "n_errors": self.n_errors,
            "n_warnings": self.n_warnings,
            "measurement_files": [item.to_dict() for item in self.measurement_files],
            "classifiers": [item.to_dict() for item in self.classifiers],
            "geojson_files": [item.to_dict() for item in self.geojson_files],
            "messages": [item.to_dict() for item in self.messages],
        }


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_report(report: CheckReport) -> None:
    report.output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = report.output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    (report.output_dir / "check_report.json").write_text(
        json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    lines = [
        "QUXYCell check report",
        f"Project: {report.project_dir}",
        f"Output: {report.output_dir}",
        f"Status: {'PASS' if report.ok else 'FAIL'}",
        f"Errors: {report.n_errors}",
        f"Warnings: {report.n_warnings}",
        "",
        f"Measurement files: {len(report.measurement_files)}",
        f"Classifier JSON files: {len(report.classifiers)}",
        f"Simple classifiers: {sum(1 for item in report.classifiers if item.is_simple)}",
        f"GeoJSON files: {len(report.geojson_files)}",
        "",
        "Messages:",
    ]
    if report.messages:
        for message in report.messages:
            suffix = f" [{message.path}]" if message.path else ""
            lines.append(f"- {message.level.upper()} {message.code}: {message.message}{suffix}")
    else:
        lines.append("- No errors or warnings.")
    (report.output_dir / "check_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    _write_csv(
        tables_dir / "measurement_files.csv",
        [
            {
                "path": str(item.path),
                "delimiter": "tab" if item.delimiter == "\t" else "comma",
                "n_columns": item.n_columns,
                "n_rows": item.n_rows if item.n_rows is not None else "",
                "columns": "|".join(item.columns),
            }
            for item in report.measurement_files
        ],
        ["path", "delimiter", "n_columns", "n_rows", "columns"],
    )
    _write_csv(
        tables_dir / "classifier_report.csv",
        [
            {
                "path": str(item.path),
                "name": item.name,
                "measurement_column": item.measurement_column or "",
                "threshold": "" if item.threshold is None else item.threshold,
                "is_simple": item.is_simple,
                "reason": item.reason,
            }
            for item in report.classifiers
        ],
        ["path", "name", "measurement_column", "threshold", "is_simple", "reason"],
    )
    _write_csv(
        tables_dir / "geojson_report.csv",
        [
            {
                "path": str(item.path),
                "readable": item.readable,
                "n_features": item.n_features if item.n_features is not None else "",
                "object_type_counts": json.dumps(item.object_type_counts, sort_keys=True),
                "class_counts": json.dumps(item.class_counts, sort_keys=True),
                "name_counts": json.dumps(item.name_counts, sort_keys=True),
                "error": item.error,
            }
            for item in report.geojson_files
        ],
        [
            "path",
            "readable",
            "n_features",
            "object_type_counts",
            "class_counts",
            "name_counts",
            "error",
        ],
    )
    _write_csv(
        tables_dir / "validation_messages.csv",
        [item.to_dict() for item in report.messages],
        ["level", "code", "message", "path"],
    )


def check(
    project_dir: str | Path,
    output_dir: str | Path = "outputs/qxy_check",
    *,
    count_rows: bool = False,
) -> CheckReport:
    """Inspect and validate a manually exported QuPath project folder.

    The function writes a report folder and returns the same information as a Python object.
    """

    project_path = Path(project_dir).expanduser().resolve()
    output_path = Path(output_dir).expanduser().resolve()
    messages: list[Message] = []

    if not project_path.exists():
        messages.append(
            Message(
                level="error",
                code="project.missing",
                message=f"Project directory does not exist: {project_path}",
                path=str(project_path),
            )
        )
        report = CheckReport(project_path, output_path, [], [], [], messages)
        _write_report(report)
        return report

    measurement_paths = discover_measurement_files(project_path)
    measurement_files: list[MeasurementFile] = []
    for path in measurement_paths:
        try:
            measurement_files.append(summarize_measurement_file(path, count_rows=count_rows))
        except Exception as exc:
            messages.append(
                Message(
                    level="error",
                    code="measurements.unreadable",
                    message=f"Measurement file could not be read: {exc}",
                    path=str(path),
                )
            )
    messages.extend(validate_measurement_files(measurement_files))

    classifier_paths = discover_classifier_files(project_path)
    classifiers = parse_classifiers(classifier_paths)
    messages.extend(validate_classifiers(classifiers, measurement_files))

    geojson_paths = discover_geojson_files(project_path)
    geojson_files = summarize_geojson_files(geojson_paths)
    messages.extend(validate_geojson_files(geojson_files))

    report = CheckReport(
        project_dir=project_path,
        output_dir=output_path,
        measurement_files=measurement_files,
        classifiers=classifiers,
        geojson_files=geojson_files,
        messages=messages,
    )
    _write_report(report)
    return report
