import csv
import json
from pathlib import Path

import qxycell as qxy


def _write_measurements(project_dir: Path) -> None:
    with (project_dir / "measurements.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Image",
                "Object ID",
                "Centroid X µm",
                "Centroid Y µm",
                "Cell: CD3: Mean",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "Image": "sample_A.tif",
                "Object ID": "1",
                "Centroid X µm": "10",
                "Centroid Y µm": "20",
                "Cell: CD3: Mean": "0.8",
            }
        )


def _write_classifier(project_dir: Path, threshold: float) -> None:
    classifier_dir = project_dir / "classifiers" / "object_classifiers"
    classifier_dir.mkdir(parents=True)
    (classifier_dir / "CD3.json").write_text(
        json.dumps(
            {
                "function": {
                    "classifier_fun": "ClassifyByMeasurementFunction",
                    "measurement": "Cell: CD3: Mean",
                    "threshold": threshold,
                }
            }
        ),
        encoding="utf-8",
    )


def _read_threshold_cell(path: Path, image: str = "sample_A.tif") -> str:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0][image]


def test_generate_threshold_table_writes_fresh_timestamped_table_from_json(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_measurements(project_dir)
    _write_classifier(project_dir, threshold=0.42)

    first = qxy.generate_threshold_table(project_dir, output_dir=output_dir)
    second = qxy.generate_threshold_table(project_dir, output_dir=output_dir)

    assert first.exists()
    assert second.exists()
    assert first != second
    assert first.name.startswith("thresholds_")
    assert second.name.startswith("thresholds_")
    assert _read_threshold_cell(first) == "0.42"


def test_check_uses_existing_threshold_table_instead_of_refreshing_from_json(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_measurements(project_dir)
    _write_classifier(project_dir, threshold=0.42)
    (project_dir / "thresholds.tsv").write_text(
        "compartment\tmarker\tmeasurement_column\tsample_A.tif\n"
        "Cell\tCD3\tCell: CD3: Mean\t0.9\n",
        encoding="utf-8",
    )

    report = qxy.check(project_dir, output_dir=output_dir)

    assert report.ok
    assert len(report.classifiers) == 1
    assert report.classifiers[0].threshold == 0.9
    assert report.classifiers[0].path.name.startswith("thresholds.tsv#")
    assert report.active_threshold_source == (project_dir / "thresholds.tsv").resolve()
    assert report.active_threshold_source_kind == "manual_threshold_file"
    assert "Threshold source: " in str(report)
    assert str((project_dir / "thresholds.tsv").resolve()) in str(report)
    report_text = report.report_path.read_text(encoding="utf-8")
    assert f"Active threshold source: {(project_dir / 'thresholds.tsv').resolve()}" in report_text
    assert "Generated threshold template:" in report_text
    generated_tables = sorted((output_dir / "tables").glob("thresholds_*.tsv"))
    assert generated_tables
    assert report.generated_threshold_template == generated_tables[-1]
    assert _read_threshold_cell(generated_tables[-1]) == "0.42"


def test_check_generates_and_uses_threshold_table_when_no_table_exists(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_measurements(project_dir)
    _write_classifier(project_dir, threshold=0.42)

    report = qxy.check(project_dir, output_dir=output_dir)

    assert report.ok
    assert len(report.classifiers) == 1
    assert report.classifiers[0].threshold == 0.42
    assert "tables" in report.classifiers[0].path.parts
    assert report.classifiers[0].path.name.startswith("thresholds_")
    assert report.active_threshold_source == report.generated_threshold_template
    assert report.active_threshold_source_kind == "generated_threshold_template"
    assert len(list((output_dir / "tables").glob("thresholds_*.tsv"))) == 1


def test_check_can_use_explicit_threshold_file(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_measurements(project_dir)
    _write_classifier(project_dir, threshold=0.42)
    (project_dir / "thresholds.tsv").write_text(
        "compartment\tmarker\tmeasurement_column\tsample_A.tif\n"
        "Cell\tCD3\tCell: CD3: Mean\t0.9\n",
        encoding="utf-8",
    )
    explicit = project_dir / "custom_thresholds.tsv"
    explicit.write_text(
        "compartment\tmarker\tmeasurement_column\tsample_A.tif\n"
        "Cell\tCD3\tCell: CD3: Mean\t0.7\n",
        encoding="utf-8",
    )

    report = qxy.check(project_dir, output_dir=output_dir, threshold_file=explicit)

    assert report.ok
    assert report.classifiers[0].threshold == 0.7
    assert report.classifiers[0].path.name.startswith("custom_thresholds.tsv#")
    assert report.active_threshold_source == explicit.resolve()
    assert report.active_threshold_source_kind == "manual_threshold_file"


def test_check_ignores_blank_manual_threshold_rows_without_warnings(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_measurements(project_dir)
    (project_dir / "thresholds.tsv").write_text(
        "compartment\tmarker\tmeasurement_column\tsample_A.tif\n"
        "Cell\tCD3\tCell: CD3: Mean\t0.9\n"
        "Cytoplasm\tCD3\tCytoplasm: CD3: Mean\t\n"
        "Membrane\tCD3\tMembrane: CD3: Mean\t#N/A\n"
        "Nucleus\tCD3\tNucleus: CD3: Mean\tn/a\n",
        encoding="utf-8",
    )

    report = qxy.check(project_dir, output_dir=output_dir)

    invalid_manual_rows = [
        message
        for message in report.messages
        if message.code == "classifiers.unsupported"
        and "invalid_manual_threshold_row" in message.message
    ]
    assert report.ok
    assert len(report.classifiers) == 1
    assert report.classifiers[0].measurement_column == "Cell: CD3: Mean"
    assert report.classifiers[0].threshold == 0.9
    assert invalid_manual_rows == []
