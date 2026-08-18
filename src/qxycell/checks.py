"""Project check/preflight workflow."""

from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Any

from qxycell.classifiers import (
    classifier_threshold_conflicts,
    compartment_from_measurement_column,
    discover_classifier_files,
    discover_threshold_files,
    marker_name_from_measurement_column,
    measurement_columns_for_threshold_template,
    parse_classifiers,
    parse_threshold_files,
    select_threshold_file,
    unresolved_threshold_conflicts,
    validate_classifiers,
)
from qxycell.geojson import discover_geojson_files, summarize_geojson_files, validate_geojson_files
from qxycell.markers import marker_name_from_classifier_name
from qxycell.measurements import (
    MEASUREMENT_TEXT_ENCODING,
    discover_measurement_files,
    summarize_measurement_file,
    validate_measurement_files,
)
from qxycell.paths import OUTPUT_TIMESTAMP_FORMAT, resolve_output_dir
from qxycell.types import ClassifierDefinition, GeoJsonFile, MeasurementFile, Message


def _safe_annotation_column(label: str) -> str:
    safe = "".join(ch if ch.isalnum() else "_" for ch in str(label).strip())
    safe = "_".join(part for part in safe.split("_") if part)
    return f"annotation__{safe or 'Unclassified'}"


@dataclass(frozen=True)
class CheckReport:
    """Combined inventory and validation report for a QuPath project folder."""

    project_dir: Path
    output_dir: Path
    measurement_files: list[MeasurementFile]
    classifiers: list[ClassifierDefinition]
    geojson_files: list[GeoJsonFile]
    messages: list[Message]
    measurement_core_counts: dict[str, int] | None = None
    geojson_core_annotation_counts: dict[str, int] | None = None
    active_threshold_source: Path | None = None
    active_threshold_source_kind: str = "none"
    generated_threshold_template: Path | None = None
    classifier_conflicts: list[dict[str, Any]] = field(default_factory=list)

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

    def geojson_object_count(self, object_type: str) -> int:
        """Return the total count of a GeoJSON objectType across discovered files."""

        target = object_type.lower()
        total = 0
        for geojson_file in self.geojson_files:
            for observed_type, count in geojson_file.object_type_counts.items():
                if observed_type.lower() == target:
                    total += int(count)
        return total

    @property
    def n_annotation_features(self) -> int:
        return self.geojson_object_count("annotation")

    @property
    def n_tma_core_features(self) -> int:
        return self.geojson_object_count("tmaCore")

    @property
    def n_cell_features(self) -> int:
        return self.geojson_object_count("cell")

    @property
    def annotation_label_counts(self) -> dict[str, int]:
        """All annotation names/classifications discovered in GeoJSON files."""
        counts: Counter[str] = Counter()
        for geojson_file in self.geojson_files:
            for object_type, labels in geojson_file.labels_by_object_type.items():
                if object_type.lower() != "annotation":
                    continue
                for label, count in labels.items():
                    if label and label.lower() not in {"", "none", "null"}:
                        counts[label] += int(count)
        return dict(sorted(counts.items()))

    @property
    def annotation_obs_columns(self) -> dict[str, list[str]]:
        """Map expected AnnData columns to the annotation labels feeding them."""
        mapping: dict[str, list[str]] = {}
        for label in self.annotation_label_counts:
            column = "Sample" if "sample" in label.lower() else _safe_annotation_column(label)
            mapping.setdefault(column, []).append(label)
        return {
            column: sorted(labels)
            for column, labels in sorted(mapping.items())
        }

    @property
    def annotation_assignments(self) -> list[dict[str, Any]]:
        """Planned source annotation to AnnData destination assignments."""
        assignments = []
        for label, count in self.annotation_label_counts.items():
            is_sample = "sample" in label.lower()
            assignments.append(
                {
                    "source_annotation": label,
                    "n_geojson_features": int(count),
                    "destination_column": "Sample" if is_sample else _safe_annotation_column(label),
                    "destination_value": label if is_sample else True,
                    "assignment_type": "categorical_value" if is_sample else "boolean_membership",
                }
            )
        return assignments

    @property
    def n_measurement_core_labels(self) -> int:
        return len(self.measurement_core_counts or {})

    @property
    def n_measurement_core_cells(self) -> int:
        return sum(int(count) for count in (self.measurement_core_counts or {}).values())

    @property
    def n_geojson_core_annotation_features(self) -> int:
        return sum(
            int(count) for count in (self.geojson_core_annotation_counts or {}).values()
        )

    @property
    def geojson_tma_core_counts(self) -> dict[str, int]:
        counts: Counter[str] = Counter(self.geojson_core_annotation_counts or {})
        for geojson_file in self.geojson_files:
            for object_type, labels in geojson_file.labels_by_object_type.items():
                if object_type.lower() != "tmacore":
                    continue
                for label, count in labels.items():
                    if label and label.lower() not in {"", "none", "null"}:
                        counts[label] += int(count)
        return dict(sorted(counts.items()))

    @property
    def n_geojson_tma_core_ids(self) -> int:
        return len(self.geojson_tma_core_counts)

    @property
    def n_geojson_tma_core_features(self) -> int:
        return sum(int(count) for count in self.geojson_tma_core_counts.values())

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

        annotation_labels = ", ".join(
            f"{label} ({count})" for label, count in self.annotation_label_counts.items()
        ) or "none"
        annotation_columns = ", ".join(self.annotation_obs_columns) or "none"

        return [
            f"QXYCell check: {'PASS' if self.ok else 'FAIL'}",
            f"Errors: {self.n_errors}",
            f"Warnings: {self.n_warnings}",
            f"Measurement files: {len(self.measurement_files)}",
            f"Classifier definitions: {len(self.classifiers)}",
            f"Simple classifiers: {sum(1 for item in self.classifiers if item.is_simple)}",
            f"Threshold source: {_display_threshold_source(self)}",
            f"Conflicting classifier channels: {len(self.classifier_conflicts)}",
            "Thresholds applied: no (check() validates definitions only)",
            "Cell typing applied: no (check() does not load cell type logic)",
            "LLM prompt generated: no",
            f"GeoJSON files: {len(self.geojson_files)}",
            f"Annotation features: {self.n_annotation_features}",
            f"Annotation names: {annotation_labels}",
            f"AnnData annotation columns: {annotation_columns}",
            f"Measurement CoreID values: {self.n_measurement_core_labels}",
            f"Cell features: {self.n_cell_features}",
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
            "active_threshold_source": (
                str(self.active_threshold_source) if self.active_threshold_source else None
            ),
            "active_threshold_source_kind": self.active_threshold_source_kind,
            "generated_threshold_template": (
                str(self.generated_threshold_template)
                if self.generated_threshold_template
                else None
            ),
            "classifier_conflicts": self.classifier_conflicts,
            "processing": {
                "thresholds_applied": False,
                "threshold_application_source": None,
                "celltyping_applied": False,
                "celltype_logic_source": None,
                "llm_prompt_generated": False,
                "llm_prompt_path": None,
            },
            "geojson_files": [item.to_dict() for item in self.geojson_files],
            "geojson_object_counts": {
                "annotation": self.n_annotation_features,
                "tmaCore": self.n_tma_core_features,
                "cell": self.n_cell_features,
            },
            "annotation_label_counts": self.annotation_label_counts,
            "annotation_obs_columns": self.annotation_obs_columns,
            "annotation_assignments": self.annotation_assignments,
            "measurement_core_counts": dict(self.measurement_core_counts or {}),
            "geojson_core_annotation_counts": dict(
                self.geojson_core_annotation_counts or {}
            ),
            "geojson_tma_core_counts": self.geojson_tma_core_counts,
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
    conflicts = classifier_threshold_conflicts(classifiers)
    conflicts_by_measurement = {
        str(item["measurement_column"]): item
        for item in conflicts
        if item.get("image") is None
    }
    marker_lookup = _marker_lookup_by_measurement_column(classifiers)
    rows = [
        {
            "compartment": compartment_from_measurement_column(column),
            "marker": marker_lookup.get(
                column,
                marker_name_from_measurement_column(column),
            ),
            "measurement_column": column,
            "classifier_conflict": column in conflicts_by_measurement,
            "candidate_classifiers": "|".join(
                conflicts_by_measurement.get(column, {}).get("classifier_names", [])
            ),
            "candidate_thresholds": "|".join(
                str(value)
                for value in conflicts_by_measurement.get(column, {}).get("thresholds", [])
            ),
            "candidate_sources": "|".join(
                conflicts_by_measurement.get(column, {}).get("sources", [])
            ),
            **{
                image: (
                    ""
                    if column in conflicts_by_measurement
                    else _prefill_threshold(threshold_lookup, image, column)
                )
                for image in images
            },
        }
        for column in measurement_columns_for_threshold_template(measurement_files)
    ]
    return rows, [
        "compartment",
        "marker",
        "measurement_column",
        "classifier_conflict",
        "candidate_classifiers",
        "candidate_thresholds",
        "candidate_sources",
        *images,
    ]


def _marker_lookup_by_measurement_column(
    classifiers: list[ClassifierDefinition],
) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for classifier in classifiers:
        if not classifier.is_simple or classifier.measurement_column is None:
            continue
        lookup.setdefault(
            str(classifier.measurement_column),
            marker_name_from_classifier_name(classifier.name),
        )
    return lookup


def _write_threshold_table(
    thresholds_dir: Path,
    measurement_files: list[MeasurementFile],
    classifiers: list[ClassifierDefinition],
    *,
    always_timestamped: bool,
    filename: str | None = None,
) -> Path:
    """Write a TSV threshold table."""

    template_path = thresholds_dir / (
        filename
        if filename is not None
        else _threshold_template_filename(
            thresholds_dir,
            always_timestamped=always_timestamped,
        )
    )
    rows, columns = _threshold_table_rows(measurement_files, classifiers)
    _write_csv(
        template_path,
        rows,
        columns,
        delimiter="\t",
    )
    conflicts = classifier_threshold_conflicts(classifiers)
    _write_csv(
        thresholds_dir.parent / "tables" / "classifier_conflicts.csv",
        [
            {
                "measurement_column": item["measurement_column"],
                "image": item.get("image") or "",
                "n_definitions": item["n_definitions"],
                "classifier_names": "|".join(item["classifier_names"]),
                "thresholds": "|".join(str(value) for value in item["thresholds"]),
                "sources": "|".join(item["sources"]),
            }
            for item in conflicts
        ],
        [
            "measurement_column",
            "image",
            "n_definitions",
            "classifier_names",
            "thresholds",
            "sources",
        ],
    )
    return template_path


def threshold_tables_dir(
    project_dir: str | Path,
    output_dir: str | Path | None = None,
) -> Path:
    """Return the generated threshold-table folder for a QXYCell output folder."""

    project_path = Path(project_dir).expanduser().resolve()
    output_path = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else resolve_output_dir(None, project_dir=project_path, project_output_kind="run")
    )
    return output_path / "thresholds"


def _threshold_template_filename(thresholds_dir: Path, *, always_timestamped: bool) -> str:
    timestamp = datetime.now().strftime(OUTPUT_TIMESTAMP_FORMAT)
    filename = f"thresholds_{timestamp}.tsv" if always_timestamped else "thresholds.tsv"
    path = thresholds_dir / filename
    if not always_timestamped or not path.exists():
        return filename
    stem = path.stem
    suffix = path.suffix
    counter = 2
    while (thresholds_dir / f"{stem}-{counter}{suffix}").exists():
        counter += 1
    return f"{stem}-{counter}{suffix}"


def _unique_measurement_images(measurement_files: list[MeasurementFile]) -> list[str]:
    images: list[str] = []
    for measurement_file in measurement_files:
        if "Image" not in measurement_file.columns:
            continue
        try:
            with measurement_file.path.open(
                newline="",
                errors="replace",
                encoding=MEASUREMENT_TEXT_ENCODING,
            ) as handle:
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


def _display_threshold_source(report: CheckReport) -> str:
    if report.active_threshold_source is not None:
        return str(report.active_threshold_source)
    if report.active_threshold_source_kind == "object_classifiers":
        return "object classifier JSON files"
    if report.active_threshold_source_kind == "none":
        return "none"
    return report.active_threshold_source_kind


def _valid_core_label(value: Any) -> str:
    label = str(value).strip()
    if label.lower() in {
        "",
        "nan",
        "none",
        "<na>",
        "na",
        "null",
        "unassigned",
        "root object (image)",
    }:
        return ""
    return label


def _measurement_core_counts(
    measurement_files: list[MeasurementFile],
    source_cols: tuple[str, ...] = ("TMA Core",),
) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for measurement_file in measurement_files:
        available_cols = [column for column in source_cols if column in measurement_file.columns]
        if not available_cols:
            continue
        try:
            with measurement_file.path.open(
                newline="",
                errors="replace",
                encoding=MEASUREMENT_TEXT_ENCODING,
            ) as handle:
                reader = csv.DictReader(handle, delimiter=measurement_file.delimiter)
                for row in reader:
                    for column in available_cols:
                        label = _valid_core_label(row.get(column, ""))
                        if label:
                            counts[label] += 1
                            break
        except Exception:
            continue
    return dict(sorted(counts.items()))


def _geojson_core_annotation_counts(
    geojson_files: list[GeoJsonFile],
    measurement_core_counts: dict[str, int],
) -> dict[str, int]:
    measurement_labels = set(measurement_core_counts)
    if not measurement_labels:
        return {}
    counts: Counter[str] = Counter()
    for geojson_file in geojson_files:
        for object_type, labels in geojson_file.labels_by_object_type.items():
            if object_type.lower() != "annotation":
                continue
            for label, count in labels.items():
                if label in measurement_labels:
                    counts[label] += int(count)
    return dict(sorted(counts.items()))


def _write_report(report: CheckReport) -> None:
    report.output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = report.output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    threshold_template_path = report.generated_threshold_template

    (report.output_dir / "check_report.json").write_text(
        json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    def _add_label_count(target: dict[str, int], label: str, count: int) -> None:
        if label and label.lower() not in ("", "none", "null"):
            target[label] = target.get(label, 0) + count

    # Aggregate labels by GeoJSON object type. QuPath GeoJSON files can include
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

    measurement_core_counts = report.measurement_core_counts or {}
    # Categorise annotation labels independently from measurement CoreID.
    ignore_labels = {k: v for k, v in annotation_labels.items() if "ignore" in k.lower()}
    sample_labels = {k: v for k, v in annotation_labels.items() if "sample" in k.lower()}
    other_labels = {
        k: v for k, v in annotation_labels.items()
        if (
            k not in ignore_labels
            and k not in sample_labels
        )
    }

    def _fmt_labels(d: dict[str, int]) -> str:
        if not d:
            return "none"
        return ", ".join(f"{k} ({v})" for k, v in sorted(d.items()))

    annotation_assignment_lines = []
    for assignment in report.annotation_assignments:
        source = assignment["source_annotation"]
        count = assignment["n_geojson_features"]
        destination = assignment["destination_column"]
        if assignment["assignment_type"] == "categorical_value":
            result = f'adata.obs["{destination}"] = "{source}"'
        else:
            result = f'adata.obs["{destination}"] = True'
        annotation_assignment_lines.append(
            f"  {source} ({count} GeoJSON feature{'s' if count != 1 else ''}) -> {result}"
        )
    if not annotation_assignment_lines:
        annotation_assignment_lines = ["  none"]

    classifier_conflict_lines = []
    for conflict in report.classifier_conflicts:
        scope = conflict.get("image") or "global"
        candidates = ", ".join(
            f"{item['classifier_name']}={item['threshold']}"
            for item in conflict["candidates"]
        )
        classifier_conflict_lines.append(
            f"    {conflict['measurement_column']} [{scope}]: {candidates}"
        )
    if not classifier_conflict_lines:
        classifier_conflict_lines = ["    none"]

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
        f"  Active threshold source: {_display_threshold_source(report)}",
        f"  Generated threshold table: {threshold_template_path or 'not generated by check()'}",
        f"  Conflicting measurement channels: {len(report.classifier_conflicts)}",
        *classifier_conflict_lines,
        "Processing performed by check():",
        "  Thresholds applied: no — definitions were inspected/validated only",
        "  Cell typing applied: no — no cell type logic was loaded",
        "  LLM prompt generated: no",
        f"GeoJSON files: {len(report.geojson_files)}",
        "GeoJSON object counts:",
        f"  Annotation features: {report.n_annotation_features}",
        f"  QuPath tmaCore objects: {report.n_tma_core_features}",
        f"  Cell features     : {report.n_cell_features}",
        "Measurement CoreID column:",
        f"  Unique CoreIDs    : {report.n_measurement_core_labels}",
        f"  Assigned cells    : {report.n_measurement_core_cells}",
        f"  CoreIDs           : {_fmt_labels(measurement_core_counts)}",
        "",
        "Annotations:",
        f"  All names: {_fmt_labels(report.annotation_label_counts)}",
        f"  Sample   : {_fmt_labels(sample_labels)}",
        f"  Ignore   : {_fmt_labels(ignore_labels)}",
        f"  Other    : {_fmt_labels(other_labels)}",
        "Planned AnnData assignments from annotations:",
        *annotation_assignment_lines,
        f"TMA cores: {_fmt_labels(tma_labels)}",
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
        tables_dir / "classifier_conflicts.csv",
        [
            {
                "measurement_column": item["measurement_column"],
                "image": item.get("image") or "",
                "n_definitions": item["n_definitions"],
                "classifier_names": "|".join(item["classifier_names"]),
                "thresholds": "|".join(str(value) for value in item["thresholds"]),
                "sources": "|".join(item["sources"]),
            }
            for item in report.classifier_conflicts
        ],
        [
            "measurement_column",
            "image",
            "n_definitions",
            "classifier_names",
            "thresholds",
            "sources",
        ],
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
        tables_dir / "coreid_report.csv",
        [
            {
                "label": label,
                "measurement_cell_count": measurement_core_counts.get(label, ""),
                "status": "measurement_coreid",
            }
            for label in sorted(measurement_core_counts)
        ],
        ["label", "measurement_cell_count", "status"],
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
    ``project_dir``, writes a threshold table under ``output_dir/thresholds/``,
    and returns the table path. Existing threshold tables are never modified or
    used as input here. If ``output_dir`` is omitted, a timestamped run output
    folder is created beside the QuPath project folder.
    """

    project_path = Path(project_dir).expanduser().resolve()
    if not project_path.exists():
        raise FileNotFoundError(f"QuPath project folder does not exist: {project_path}")
    thresholds_dir = threshold_tables_dir(project_path, output_dir=output_dir)

    measurement_files = [
        summarize_measurement_file(path, count_rows=count_rows)
        for path in discover_measurement_files(project_path)
    ]
    classifiers = parse_classifiers(discover_classifier_files(project_path))
    return _write_threshold_table(
        thresholds_dir,
        measurement_files,
        classifiers,
        always_timestamped=True,
    )


def write_classifier_threshold_table(
    project_dir: str | Path,
    output_dir: str | Path,
    classifiers: list[ClassifierDefinition],
) -> Path:
    """Replace the stable table containing classifier-derived thresholds."""

    project_path = Path(project_dir).expanduser().resolve()
    measurement_files = [
        summarize_measurement_file(path, count_rows=False)
        for path in discover_measurement_files(project_path)
    ]
    return _write_threshold_table(
        threshold_tables_dir(project_path, output_dir=output_dir),
        measurement_files,
        classifiers,
        always_timestamped=False,
        filename="classifier_thresholds.tsv",
    )


def inspect_project(
    project_dir: str | Path,
    *,
    count_rows: bool = False,
    threshold_file: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> CheckReport:
    """Inspect a QuPath project folder without writing reports."""

    project_path = Path(project_dir).expanduser().resolve()
    output_path = Path(output_dir).expanduser().resolve() if output_dir is not None else project_path
    messages: list[Message] = []

    if not project_path.exists():
        messages.append(
            Message(
                level="error",
                code="project.missing",
                message=f"QuPath project folder does not exist: {project_path}",
                path=str(project_path),
            )
        )
        return CheckReport(project_path, output_path, [], [], [], messages)

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
    classifier_conflicts = classifier_threshold_conflicts(object_classifiers)
    for conflict in classifier_conflicts:
        scope = conflict.get("image") or "global"
        candidates = ", ".join(
            f"{item['classifier_name']}={item['threshold']}"
            for item in conflict["candidates"]
        )
        messages.append(
            Message(
                level="warning",
                code="classifiers.conflicting_thresholds",
                message=(
                    "Multiple thresholds target the same measurement and scope: "
                    f"{conflict['measurement_column']} [{scope}] ({candidates}). "
                    "Generated threshold tables leave this channel blank for review."
                ),
            )
        )
    generated_threshold_path = None
    active_threshold_source = None
    active_threshold_source_kind = "none"
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
        threshold_paths = discover_threshold_files(project_path, output_dir=output_path)
        threshold_path, ignored_threshold_paths = select_threshold_file(threshold_paths)

    if threshold_path is not None and threshold_path.exists():
        classifiers = parse_threshold_files([threshold_path])
        active_threshold_source = threshold_path
        active_threshold_source_kind = "manual_threshold_file"
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
        active_threshold_source = generated_threshold_path
        active_threshold_source_kind = "generated_threshold_template"
    else:
        classifiers = object_classifiers
        if classifiers:
            active_threshold_source_kind = "object_classifiers"
    classifier_messages = validate_classifiers(classifiers, measurement_files)
    unresolved_conflicts = (
        unresolved_threshold_conflicts(threshold_path)
        if threshold_path is not None and threshold_path.exists()
        else []
    )
    if unresolved_conflicts:
        classifier_messages = [
            message
            for message in classifier_messages
            if message.code not in {
                "classifiers.no_simple_thresholds",
                "classifiers.missing",
            }
        ]
        messages.append(
            Message(
                level="warning",
                code="classifiers.unresolved_conflicts",
                message=(
                    f"{len(unresolved_conflicts)} conflict-marked threshold row(s) "
                    "still require per-image values before thresholding can be applied."
                ),
                path=str(threshold_path),
            )
        )
    messages.extend(classifier_messages)

    geojson_paths = discover_geojson_files(project_path)
    geojson_files = summarize_geojson_files(geojson_paths)
    messages.extend(validate_geojson_files(geojson_files))
    measurement_core_counts = _measurement_core_counts(measurement_files)
    geojson_core_annotation_counts = _geojson_core_annotation_counts(
        geojson_files,
        measurement_core_counts,
    )

    report = CheckReport(
        project_dir=project_path,
        output_dir=output_path,
        measurement_files=measurement_files,
        classifiers=classifiers,
        geojson_files=geojson_files,
        messages=messages,
        measurement_core_counts=measurement_core_counts,
        geojson_core_annotation_counts=geojson_core_annotation_counts,
        active_threshold_source=active_threshold_source,
        active_threshold_source_kind=active_threshold_source_kind,
        generated_threshold_template=generated_threshold_path,
        classifier_conflicts=classifier_conflicts,
    )
    return report


def check(
    project_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    count_rows: bool = False,
    threshold_file: str | Path | None = None,
) -> CheckReport:
    """Inspect and validate a QuPath project folder.

    The function writes a report folder and returns the same information as a Python object.
    """

    project_path = Path(project_dir).expanduser().resolve()
    output_path = resolve_output_dir(
        output_dir,
        project_dir=project_path,
        project_output_kind="check",
    )
    report = inspect_project(
        project_path,
        count_rows=count_rows,
        threshold_file=threshold_file,
        output_dir=output_path,
    )
    _write_report(report)
    return report
