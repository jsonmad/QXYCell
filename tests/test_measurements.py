import json

import pytest

import qxycell as qxy
from qxycell.measurements import (
    REQUIRED_MEASUREMENT_COLUMNS,
    summarize_measurement_file,
    validate_measurement_files,
)


def _write_measurements(path, *, encoding="utf-8", columns=REQUIRED_MEASUREMENT_COLUMNS):
    delimiter = "\t" if path.suffix == ".tsv" else ","
    header = delimiter.join(columns)
    row = delimiter.join(["image_1", "cell_1", "10.0", "20.0"][: len(columns)])
    path.write_text(f"{header}\n{row}\n", encoding=encoding)


def _write_classifier(project_dir):
    classifier_dir = project_dir / "classifiers" / "object_classifiers"
    classifier_dir.mkdir(parents=True)
    (classifier_dir / "CD3.json").write_text(
        json.dumps(
            {
                "function": {
                    "classifier_fun": "ClassifyByMeasurementFunction",
                    "measurement": "Centroid X µm",
                    "threshold": 0.5,
                }
            }
        ),
        encoding="utf-8",
    )


@pytest.mark.parametrize("suffix", [".tsv", ".csv"])
def test_utf8_measurement_headers_validate(suffix, tmp_path):
    path = tmp_path / f"measurements{suffix}"
    _write_measurements(path, encoding="utf-8")

    summary = summarize_measurement_file(path)

    assert summary.columns == REQUIRED_MEASUREMENT_COLUMNS
    assert validate_measurement_files([summary]) == []


def test_utf8_bom_measurement_header_validates_without_bom_in_image(tmp_path):
    path = tmp_path / "measurements.tsv"
    _write_measurements(path, encoding="utf-8-sig")

    summary = summarize_measurement_file(path)

    assert summary.columns[0] == "Image"
    assert not summary.columns[0].startswith("\ufeff")
    assert validate_measurement_files([summary]) == []


def test_missing_measurement_column_validation_still_fails(tmp_path):
    path = tmp_path / "measurements.tsv"
    columns = tuple(
        column for column in REQUIRED_MEASUREMENT_COLUMNS if column != "Centroid Y µm"
    )
    _write_measurements(path, columns=columns)

    messages = validate_measurement_files([summarize_measurement_file(path)])

    assert len(messages) == 1
    assert messages[0].code == "measurements.required_columns_missing"
    assert "Centroid Y µm" in messages[0].message


@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig"])
def test_check_accepts_utf8_measurement_headers(encoding, tmp_path):
    project_dir = tmp_path / "qupath_project"
    project_dir.mkdir()
    _write_measurements(project_dir / "measurements.tsv", encoding=encoding)
    _write_classifier(project_dir)

    report = qxy.check(
        project_dir,
        output_dir=tmp_path / "check_output",
        count_rows=True,
    )

    measurement_errors = [
        message
        for message in report.messages
        if message.code == "measurements.required_columns_missing"
    ]
    assert measurement_errors == []
    assert report.measurement_files[0].columns[0] == "Image"
