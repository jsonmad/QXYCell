import csv
import json
from pathlib import Path

import qxycell as qxy
import pandas as pd
import pytest


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


def _write_conflicting_classifier(project_dir: Path, name: str, threshold: float) -> None:
    classifier_dir = project_dir / "classifiers" / "object_classifiers"
    classifier_dir.mkdir(parents=True, exist_ok=True)
    (classifier_dir / f"{name}.json").write_text(
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


def _write_two_image_measurements(project_dir: Path) -> None:
    pd.DataFrame(
        {
            "Image": ["image_A.tif", "image_A.tif", "image_B.tif", "image_B.tif"],
            "Object ID": ["1", "2", "3", "4"],
            "Centroid X µm": [0.0, 1.0, 0.0, 1.0],
            "Centroid Y µm": [0.0, 1.0, 0.0, 1.0],
            "Cell: CD3: Mean": [5.0, 15.0, 5.0, 15.0],
        }
    ).to_csv(project_dir / "measurements.csv", index=False)


def _write_alias_measurements(project_dir: Path) -> None:
    with (project_dir / "measurements.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Image",
                "Object ID",
                "Centroid X µm",
                "Centroid Y µm",
                "Cell: #945;SMA - TRITC: Mean",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "Image": "sample_A.tif",
                "Object ID": "1",
                "Centroid X µm": "10",
                "Centroid Y µm": "20",
                "Cell: #945;SMA - TRITC: Mean": "0.8",
            }
        )


def _write_alias_classifier(project_dir: Path, threshold: float) -> None:
    classifier_dir = project_dir / "classifiers" / "object_classifiers"
    classifier_dir.mkdir(parents=True)
    (classifier_dir / "aSMA.json").write_text(
        json.dumps(
            {
                "function": {
                    "classifier_fun": "ClassifyByMeasurementFunction",
                    "measurement": "Cell: #945;SMA - TRITC: Mean",
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
    assert first.parent == output_dir / "thresholds"
    assert second.parent == output_dir / "thresholds"
    assert first.name.startswith("thresholds_")
    assert second.name.startswith("thresholds_")
    assert _read_threshold_cell(first) == "0.42"


def test_generate_threshold_table_uses_classifier_filename_for_marker_name(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_alias_measurements(project_dir)
    _write_alias_classifier(project_dir, threshold=0.42)

    table_path = qxy.generate_threshold_table(project_dir, output_dir=output_dir)

    assert table_path.parent == output_dir / "thresholds"
    with table_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["marker"] == "aSMA"
    assert rows[0]["measurement_column"] == "Cell: #945;SMA - TRITC: Mean"
    assert rows[0]["sample_A.tif"] == "0.42"


def test_conflicting_classifiers_are_reported_and_generated_cells_are_blank(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_measurements(project_dir)
    _write_conflicting_classifier(project_dir, "CD3", 0.42)
    _write_conflicting_classifier(project_dir, "CD3_2", 0.90)

    report = qxy.check(project_dir, output_dir=output_dir)
    table_path = qxy.generate_threshold_table(project_dir, output_dir=output_dir)
    table = pd.read_csv(table_path, sep="\t", dtype=str).fillna("")
    conflict_table = pd.read_csv(output_dir / "tables" / "classifier_conflicts.csv")

    assert len(report.classifier_conflicts) == 1
    assert report.classifier_conflicts[0]["thresholds"] == [0.42, 0.9]
    assert any(message.code == "classifiers.conflicting_thresholds" for message in report.messages)
    assert "CD3=0.42" in report.report_path.read_text(encoding="utf-8")
    assert table.loc[0, "classifier_conflict"] == "True"
    assert table.loc[0, "candidate_classifiers"] == "CD3|CD3_2"
    assert table.loc[0, "candidate_thresholds"] == "0.42|0.9"
    assert table.loc[0, "sample_A.tif"] == ""
    assert conflict_table.loc[0, "measurement_column"] == "Cell: CD3: Mean"


def test_per_image_conflict_table_can_be_completed_and_reapplied_with_celltyping(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_two_image_measurements(project_dir)
    _write_conflicting_classifier(project_dir, "CD3", 4.0)
    _write_conflicting_classifier(project_dir, "CD3_2", 12.0)

    adata = qxy.run(project_dir, output_dir=output_dir, verbose=False)
    table_path = next((output_dir / "thresholds").glob("thresholds_*.tsv"))

    with pytest.raises(ValueError, match="unresolved classifier conflicts"):
        qxy.threshold(adata, threshold_file=table_path, verbose=False)

    table = pd.read_csv(table_path, sep="\t", dtype=str).fillna("")
    table.loc[table["measurement_column"] == "Cell: CD3: Mean", "image_A.tif"] = "10"
    table.loc[table["measurement_column"] == "Cell: CD3: Mean", "image_B.tif"] = "20"
    table.to_csv(table_path, sep="\t", index=False)

    summary = qxy.threshold(adata, threshold_file=table_path, verbose=False)
    assert summary["pos_columns"] == ["CD3_pos"]
    assert adata.obs["CD3_pos"].tolist() == [0, 1, 0, 0]

    logic = {"rules": [{"name": "CD3+", "positive": ["CD3"]}]}
    qxy.celltype(adata, logic, verbose=False)
    assert adata.obs["celltype"].tolist() == ["Unknown", "CD3+", "Unknown", "Unknown"]

    table.loc[table["measurement_column"] == "Cell: CD3: Mean", "image_A.tif"] = "0"
    table.loc[table["measurement_column"] == "Cell: CD3: Mean", "image_B.tif"] = "10"
    table.to_csv(table_path, sep="\t", index=False)
    qxy.threshold(adata, threshold_file=table_path, verbose=False)
    assert "celltype" not in adata.obs.columns
    assert "celltype__stale_celltype" in adata.obs.columns
    assert adata.obs["CD3_pos"].tolist() == [1, 1, 0, 1]

    qxy.celltype(adata, logic, verbose=False)
    assert adata.obs["celltype"].tolist() == ["CD3+", "CD3+", "Unknown", "CD3+"]


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
    assert "Generated threshold table: not generated by check()" in report_text
    assert report.generated_threshold_template is None
    assert not (output_dir / "thresholds").exists()


def test_check_uses_object_classifiers_when_no_threshold_table_exists(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_measurements(project_dir)
    _write_classifier(project_dir, threshold=0.42)

    report = qxy.check(project_dir, output_dir=output_dir)

    assert report.ok
    assert len(report.classifiers) == 1
    assert report.classifiers[0].threshold == 0.42
    assert report.classifiers[0].path.name == "CD3.json"
    assert report.active_threshold_source is None
    assert report.generated_threshold_template is None
    assert report.active_threshold_source_kind == "object_classifiers"
    assert not (output_dir / "thresholds").exists()


def test_check_does_not_reuse_threshold_table_from_old_explicit_output(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = project_dir / "outputs_260713-141927"
    stale_dir = project_dir / "outputs_260713-130128" / "thresholds"
    project_dir.mkdir()
    stale_dir.mkdir(parents=True)
    _write_measurements(project_dir)
    _write_classifier(project_dir, threshold=0.42)
    (stale_dir / "thresholds_260713-1301.tsv").write_text(
        "marker\tmeasurement_column\tthreshold\n"
        "Marker\tMarker: Mean\t0.99\n",
        encoding="utf-8",
    )

    report = qxy.check(project_dir, output_dir=output_dir)

    assert report.active_threshold_source is None
    assert report.active_threshold_source_kind == "object_classifiers"
    assert report.classifiers[0].threshold == 0.42


def test_check_uses_output_thresholds_but_ignores_generated_inventory_tables(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = project_dir / "outputs"
    tables_dir = output_dir / "tables"
    thresholds_dir = output_dir / "thresholds"
    tables_dir.mkdir(parents=True)
    thresholds_dir.mkdir(parents=True)
    _write_measurements(project_dir)
    _write_classifier(project_dir, threshold=0.42)
    (tables_dir / "measurement_files.csv").write_text(
        "path,n_columns\nmeasurements.csv,5\n",
        encoding="utf-8",
    )
    (thresholds_dir / "thresholds_250101-0101.tsv").write_text(
        "compartment\tmarker\tmeasurement_column\tsample_A.tif\n"
        "Cell\tCD3\tCell: CD3: Mean\t0.99\n",
        encoding="utf-8",
    )

    report = qxy.check(project_dir, output_dir=output_dir)

    assert report.ok
    assert [item.path.name for item in report.measurement_files] == ["measurements.csv"]
    assert report.classifiers[0].threshold == 0.99
    assert report.active_threshold_source == (thresholds_dir / "thresholds_250101-0101.tsv").resolve()
    assert report.active_threshold_source_kind == "manual_threshold_file"
    assert report.generated_threshold_template is None
    assert not any(
        message.path == str(tables_dir / "measurement_files.csv")
        for message in report.messages
    )


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


def test_threshold_archives_stale_celltype_outputs(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    _write_measurements(project_dir)
    _write_classifier(project_dir, threshold=0.42)

    adata = qxy.run(project_dir, output_dir=output_dir, verbose=False)
    adata.obs["celltype"] = "Old type"
    adata.obs["custom_celltype"] = "Old custom type"
    adata.obs["old_feature"] = 1
    adata.obs["old_derived_feature"] = 1
    adata.uns["qxycell_celltyping"] = {
        "celltype_column": "custom_celltype",
        "feature_columns": ["old_feature"],
        "derived_feature_columns": ["old_derived_feature"],
    }
    adata.uns["qxycell"]["celltyping_applied"] = True

    summary = qxy.threshold(
        adata,
        project_dir=project_dir,
        output_dir=output_dir,
        verbose=False,
    )

    assert "CD3_pos" in adata.obs.columns
    assert "celltype" not in adata.obs.columns
    assert "custom_celltype" not in adata.obs.columns
    assert "old_feature" not in adata.obs.columns
    assert "old_derived_feature" not in adata.obs.columns
    assert adata.obs["celltype__stale_celltype"].tolist() == ["Old type"]
    assert adata.obs["custom_celltype__stale_celltype"].tolist() == ["Old custom type"]
    assert adata.obs["old_feature__stale_celltype"].tolist() == [1]
    assert adata.obs["old_derived_feature__stale_celltype"].tolist() == [1]
    assert "qxycell_celltyping" not in adata.uns
    assert adata.uns["qxycell_stale_celltyping"][0]["reason"] == "thresholds_reapplied"
    assert adata.uns["qxycell_stale_celltyping"][0]["columns"] == {
        "celltype": "celltype__stale_celltype",
        "custom_celltype": "custom_celltype__stale_celltype",
        "old_feature": "old_feature__stale_celltype",
        "old_derived_feature": "old_derived_feature__stale_celltype",
    }
    assert adata.uns["qxycell"]["celltyping_applied"] is False
    assert summary["stale_celltype_columns"] == {
        "celltype": "celltype__stale_celltype",
        "custom_celltype": "custom_celltype__stale_celltype",
        "old_feature": "old_feature__stale_celltype",
        "old_derived_feature": "old_derived_feature__stale_celltype",
    }


def test_rethreshold_removes_only_stale_qxycell_owned_pos_columns(tmp_path):
    project_dir = tmp_path / "project"
    output_dir = tmp_path / "output"
    project_dir.mkdir()
    pd.DataFrame(
        {
            "Image": ["sample_A.tif"],
            "Object ID": ["1"],
            "Centroid X µm": [10.0],
            "Centroid Y µm": [20.0],
            "Cell: CD3: Mean": [0.8],
            "Cell: CD8: Mean": [0.3],
        }
    ).to_csv(project_dir / "measurements.csv", index=False)

    classifier_dir = project_dir / "classifiers" / "object_classifiers"
    classifier_dir.mkdir(parents=True)
    for marker, threshold in (("CD3", 0.5), ("CD8", 0.4)):
        (classifier_dir / f"{marker}.json").write_text(
            json.dumps(
                {
                    "function": {
                        "classifier_fun": "ClassifyByMeasurementFunction",
                        "measurement": f"Cell: {marker}: Mean",
                        "threshold": threshold,
                    }
                }
            ),
            encoding="utf-8",
        )

    initial_thresholds = tmp_path / "initial_thresholds.tsv"
    replacement_thresholds = tmp_path / "replacement_thresholds.tsv"
    pd.DataFrame(
        {
            "marker": ["CD3", "CD8"],
            "measurement_column": ["Cell: CD3: Mean", "Cell: CD8: Mean"],
            "threshold": [0.5, 0.4],
        }
    ).to_csv(initial_thresholds, sep="\t", index=False)
    pd.DataFrame(
        {
            "marker": ["CD3"],
            "measurement_column": ["Cell: CD3: Mean"],
            "threshold": [0.7],
        }
    ).to_csv(replacement_thresholds, sep="\t", index=False)

    adata = qxy.run(project_dir, output_dir=output_dir, verbose=False)
    qxy.threshold(adata, threshold_file=initial_thresholds, verbose=False)
    assert adata.uns["qxycell_thresholding"]["pos_columns"] == ["CD3_pos", "CD8_pos"]
    adata.obs["manual_pos"] = True

    summary = qxy.threshold(adata, threshold_file=replacement_thresholds, verbose=False)

    assert summary["pos_columns"] == ["CD3_pos"]
    assert "CD8_pos" not in adata.obs.columns
    assert adata.obs["manual_pos"].tolist() == [True]
