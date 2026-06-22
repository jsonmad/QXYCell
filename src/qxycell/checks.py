"""Project check/preflight workflow."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Any

from qxycell.classifiers import (
    compartment_from_measurement_column,
    discover_classifier_files,
    discover_threshold_files,
    marker_name_from_measurement_column,
    measurement_columns_for_threshold_template,
    parse_classifiers,
    parse_threshold_files,
    select_threshold_file,
    validate_classifiers,
)
from qxycell.geojson import discover_geojson_files, summarize_geojson_files, validate_geojson_files
from qxycell.measurements import (
    discover_measurement_files,
    summarize_measurement_file,
    validate_measurement_files,
)
from qxycell.paths import OUTPUT_TIMESTAMP_FORMAT, resolve_output_dir
from qxycell.types import ClassifierDefinition, GeoJsonFile, MeasurementFile, Message


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

    @property
    def report_path(self) -> Path:
        """Path to the plain-text check report written by :func:`check`."""

        return self.output_dir / "check_report.txt"

    @property
    def json_path(self) -> Path:
        """Path to the JSON check report written by :func:`check`."""

        return self.output_dir / "check_report.json"

    def summary_lines(self) -> list[str]:
        """Return a compact, notebook-friendly report summary."""

        return [
            f"QXYCell check: {'PASS' if self.ok else 'FAIL'}",
            f"Errors: {self.n_errors}",
            f"Warnings: {self.n_warnings}",
            f"Measurement files: {len(self.measurement_files)}",
            f"Classifier definitions: {len(self.classifiers)}",
            f"Simple classifiers: {sum(1 for item in self.classifiers if item.is_simple)}",
            f"GeoJSON files: {len(self.geojson_files)}",
            f"Report: {self.report_path}",
        ]

    def __str__(self) -> str:
        return "\n".join(self.summary_lines())

    def __repr__(self) -> str:
        return str(self)

    def _repr_pretty_(self, printer: Any, cycle: bool) -> None:
        if cycle:
            printer.text("CheckReport(...)")
        else:
            printer.text(str(self))

    def _repr_html_(self) -> str:
        return "<pre>" + escape(str(self)) + "</pre>"

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


def _write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    columns: list[str],
    *,
    delimiter: str = ",",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore",
            delimiter=delimiter,
        )
        writer.writeheader()
        writer.writerows(rows)


def _threshold_table_rows(
    measurement_files: list[MeasurementFile],
    classifiers: list[ClassifierDefinition],
) -> tuple[list[dict[str, Any]], list[str]]:
    images = _unique_measurement_images(measurement_files)
    threshold_lookup = _threshold_lookup(classifiers)
    rows = [
        {
            "compartment": compartment_from_measurement_column(column),
            "marker": marker_name_from_measurement_column(column),
            "measurement_column": column,
            **{
                image: _prefill_threshold(threshold_lookup, image, column)
                for image in images
            },
        }
        for column in measurement_columns_for_threshold_template(measurement_files)
    ]
    return rows, ["compartment", "marker", "measurement_column", *images]


def _write_threshold_table(
    output_dir: Path,
    measurement_files: list[MeasurementFile],
    classifiers: list[ClassifierDefinition],
    *,
    always_timestamped: bool,
) -> Path:
    """Write a TSV threshold table."""

    template_path = output_dir / "tables" / _threshold_template_filename(
        output_dir,
        always_timestamped=always_timestamped,
    )
    rows, columns = _threshold_table_rows(measurement_files, classifiers)
    _write_csv(
        template_path,
        rows,
        columns,
        delimiter="\t",
    )
    return template_path


def _write_manual_threshold_template(report: CheckReport) -> Path:
    """Write a fill-in TSV for manual marker thresholds."""

    return _write_threshold_table(
        report.output_dir,
        report.measurement_files,
        report.classifiers,
        always_timestamped=True,
    )


def _threshold_template_filename(output_dir: Path, *, always_timestamped: bool) -> str:
    timestamp = output_dir.name.removeprefix("qxy_outputs_")
    if timestamp == output_dir.name or not timestamp:
        timestamp = datetime.now().strftime(OUTPUT_TIMESTAMP_FORMAT)
    filename = f"thresholds_{timestamp}.tsv" if always_timestamped else "thresholds.tsv"
    path = output_dir / "tables" / filename
    if not always_timestamped or not path.exists():
        return filename
    stem = path.stem
    suffix = path.suffix
    counter = 2
    while (output_dir / "tables" / f"{stem}-{counter}{suffix}").exists():
        counter += 1
    return f"{stem}-{counter}{suffix}"


def _unique_measurement_images(measurement_files: list[MeasurementFile]) -> list[str]:
    images: list[str] = []
    for measurement_file in measurement_files:
        if "Image" not in measurement_file.columns:
            continue
        try:
            with measurement_file.path.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter=measurement_file.delimiter)
                for row in reader:
                    image = str(row.get("Image", "")).strip()
                    if image:
                        images.append(image)
        except Exception:
            continue
    return sorted(dict.fromkeys(images))


def _threshold_lookup(
    classifiers: list[ClassifierDefinition],
) -> dict[tuple[str | None, str], float]:
    lookup: dict[tuple[str | None, str], float] = {}
    for classifier in classifiers:
        if (
            not classifier.is_simple
            or classifier.measurement_column is None
            or classifier.threshold is None
        ):
            continue
        image = str(classifier.image).strip() if classifier.image else None
        lookup[(image, classifier.measurement_column)] = classifier.threshold
    return lookup


def _prefill_threshold(
    lookup: dict[tuple[str | None, str], float],
    image: str,
    measurement_column: str,
) -> str:
    threshold = lookup.get((image, measurement_column), lookup.get((None, measurement_column)))
    return "" if threshold is None else str(threshold)


def _output_threshold_source(report: CheckReport) -> Path | None:
    tables_dir = report.output_dir / "tables"
    for classifier in report.classifiers:
        source = Path(str(classifier.path).split("#", 1)[0])
        if source.exists() and source.parent == tables_dir:
            return source
    return None


def _latest_output_threshold_table(output_dir: Path) -> Path | None:
    tables_dir = output_dir / "tables"
    candidates = [
        path
        for path in tables_dir.glob("thresholds_*.tsv")
        if path.is_file()
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda path: (path.stat().st_mtime, str(path)))


def _write_report(report: CheckReport) -> None:
    report.output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = report.output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    threshold_template_path = _output_threshold_source(report) or _latest_output_threshold_table(
        report.output_dir
    )
    if threshold_template_path is None:
        threshold_template_path = _write_manual_threshold_template(report)

    (report.output_dir / "check_report.json").write_text(
        json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    def _add_label_count(target: dict[str, int], label: str, count: int) -> None:
        if label and label.lower() not in ("", "none", "null"):
            target[label] = target.get(label, 0) + count

    # Aggregate labels by GeoJSON object type. QuPath exports can include
    # annotation polygons, TMA cores, and individual cell objects in the same
    # file; the check report should not mix those categories.
    annotation_labels: dict[str, int] = {}
    tma_labels: dict[str, int] = {}
    cell_labels: dict[str, int] = {}
    for gf in report.geojson_files:
        for object_type, labels in gf.labels_by_object_type.items():
            object_type_lower = object_type.lower()
            for label, count in labels.items():
                if object_type_lower == "annotation":
                    _add_label_count(annotation_labels, label, count)
                elif object_type_lower == "tmacore":
                    _add_label_count(tma_labels, label, count)
                elif object_type_lower == "cell":
                    _add_label_count(cell_labels, label, count)

    # Categorise labels into known QXYCell roles.
    ignore_labels = {k: v for k, v in annotation_labels.items() if "ignore" in k.lower()}
    sample_labels = {k: v for k, v in annotation_labels.items() if "sample" in k.lower()}
    other_labels = {
        k: v for k, v in annotation_labels.items()
        if k not in ignore_labels and k not in sample_labels and k not in tma_labels
    }

    def _fmt_labels(d: dict[str, int]) -> str:
        if not d:
            return "none"
        return ", ".join(f"{k} ({v})" for k, v in sorted(d.items()))

    lines = [
        "QXYCell check report",
        f"Project: {report.project_dir}",
        f"Output: {report.output_dir}",
        f"Status: {'PASS' if report.ok else 'FAIL'}",
        f"Errors: {report.n_errors}",
        f"Warnings: {report.n_warnings}",
        "",
        f"Measurement files: {len(report.measurement_files)}",
        f"Classifier definitions: {len(report.classifiers)}",
        f"  Simple classifiers: {sum(1 for item in report.classifiers if item.is_simple)}",
        f"  Markers: {', '.join(c.name for c in report.classifiers if c.is_simple) or 'none'}",
        f"  Manual threshold template: {threshold_template_path}",
        f"GeoJSON files: {len(report.geojson_files)}",
        "",
        "Annotations:",
        f"  Sample   : {_fmt_labels(sample_labels)}",
        f"  Ignore   : {_fmt_labels(ignore_labels)}",
        f"  TMA      : {_fmt_labels(tma_labels)}",
        f"  Other    : {_fmt_labels(other_labels)}",
        f"  Cell labels: {_fmt_labels(cell_labels)}",
        "",
        "Messages:",
    ]
    if report.messages:
        for message in report.messages:
            suffix = f" [{message.path}]" if message.path else ""
            lines.append(f"- {message.level.upper()} {message.code}: {message.message}{suffix}")
    else:
        lines.append("- No errors or warnings.")
    report_text = "\n".join(lines) + "\n"
    (report.output_dir / "check_report.txt").write_text(report_text, encoding="utf-8")
    print(report_text, end="")

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
                "image": item.image or "",
                "measurement_column": item.measurement_column or "",
                "threshold": "" if item.threshold is None else item.threshold,
                "is_simple": item.is_simple,
                "reason": item.reason,
            }
            for item in report.classifiers
        ],
        ["path", "name", "image", "measurement_column", "threshold", "is_simple", "reason"],
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
                "labels_by_object_type": json.dumps(
                    item.labels_by_object_type, sort_keys=True
                ),
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
            "labels_by_object_type",
            "error",
        ],
    )
    _write_csv(
        tables_dir / "validation_messages.csv",
        [item.to_dict() for item in report.messages],
        ["level", "code", "message", "path"],
    )


def generate_threshold_table(
    project_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    count_rows: bool = False,
) -> Path:
    """Create a fresh timestamped threshold table from object classifier JSONs.

    This function reads measurement columns and object classifier JSON files from
    ``project_dir``, writes a threshold table under ``output_dir/tables/``, and
    returns the table path. Existing threshold tables are never modified or used
    as input here.
    """

    project_path = Path(project_dir).expanduser().resolve()
    if not project_path.exists():
        raise FileNotFoundError(f"Project directory does not exist: {project_path}")
    output_path = resolve_output_dir(output_dir)

    measurement_files = [
        summarize_measurement_file(path, count_rows=count_rows)
        for path in discover_measurement_files(project_path)
    ]
    classifiers = parse_classifiers(discover_classifier_files(project_path))
    return _write_threshold_table(
        output_path,
        measurement_files,
        classifiers,
        always_timestamped=True,
    )


def check(
    project_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    count_rows: bool = False,
    threshold_file: str | Path | None = None,
) -> CheckReport:
    """Inspect and validate a manually exported QuPath project folder.

    The function writes a report folder and returns the same information as a Python object.
    """

    project_path = Path(project_dir).expanduser().resolve()
    output_path = resolve_output_dir(output_dir)
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

    object_classifiers = parse_classifiers(discover_classifier_files(project_path))
    generated_threshold_path = None
    if any(classifier.is_simple for classifier in object_classifiers):
        generated_threshold_path = _write_threshold_table(
            output_path,
            measurement_files,
            object_classifiers,
            always_timestamped=True,
        )

    ignored_threshold_paths: list[Path] = []
    if threshold_file is not None:
        threshold_path = Path(threshold_file).expanduser().resolve()
        if not threshold_path.exists():
            messages.append(
                Message(
                    level="error",
                    code="classifiers.threshold_file_missing",
                    message=f"Specified threshold file does not exist: {threshold_path}",
                    path=str(threshold_path),
                )
            )
    else:
        threshold_paths = discover_threshold_files(project_path)
        threshold_path, ignored_threshold_paths = select_threshold_file(threshold_paths)

    if threshold_path is not None and threshold_path.exists():
        classifiers = parse_threshold_files([threshold_path])
        for ignored_path in ignored_threshold_paths:
            messages.append(
                Message(
                    level="warning",
                    code="classifiers.threshold_file_ignored",
                    message=(
                        "Ignoring manual threshold file because another threshold "
                        f"file was selected: {threshold_path}"
                    ),
                    path=str(ignored_path),
                )
            )
    elif generated_threshold_path is not None:
        classifiers = parse_threshold_files([generated_threshold_path])
    else:
        classifiers = object_classifiers
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
