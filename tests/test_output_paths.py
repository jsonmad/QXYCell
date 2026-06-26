import json
from pathlib import Path
import re

import pandas as pd

import qxycell as qxy
from qxycell.discovery import is_qxy_output_artifact
from qxycell.paths import latest_timestamped_output_dir, output_dir_from_adata


def _write_minimal_project(project_dir):
    project_dir.mkdir()
    pd.DataFrame(
        {
            "Image": ["sample_A.tif"],
            "Object ID": ["cell_0"],
            "Centroid X µm": [10.0],
            "Centroid Y µm": [20.0],
            "Cell: CD3: Mean": [0.8],
        }
    ).to_csv(project_dir / "measurements.tsv", sep="\t", index=False)
    classifier_dir = project_dir / "classifiers" / "object_classifiers"
    classifier_dir.mkdir(parents=True)
    (classifier_dir / "CD3.json").write_text(
        json.dumps(
            {
                "function": {
                    "classifier_fun": "ClassifyByMeasurementFunction",
                    "measurement": "Cell: CD3: Mean",
                    "threshold": 0.5,
                }
            }
        ),
        encoding="utf-8",
    )


def test_check_defaults_to_project_sibling_timestamped_output(tmp_path):
    project_dir = tmp_path / "qupath_dir"
    _write_minimal_project(project_dir)

    report = qxy.check(project_dir)

    assert report.output_dir.parent == project_dir.parent
    assert re.fullmatch(r"qupath_dir_check_\d{6}_\d{4}", report.output_dir.name)
    assert report.generated_threshold_template is None
    assert not (report.output_dir / "thresholds").exists()


def test_run_defaults_to_project_sibling_timestamped_output(tmp_path):
    project_dir = tmp_path / "qupath_dir"
    _write_minimal_project(project_dir)

    adata = qxy.run(project_dir, verbose=False)

    output_dir = output_dir_from_adata(adata)
    assert output_dir.parent == project_dir.parent
    assert re.fullmatch(r"qupath_dir_run_\d{6}_\d{4}", output_dir.name)
    assert (output_dir / "run" / "h5ad" / "qxycell.h5ad").exists()
    assert re.fullmatch(
        r"qupath_dir_check_\d{6}_\d{4}",
        Path(adata.uns["qxycell"]["check_output_dir"]).name,
    )
    assert not (output_dir / "check_report.txt").exists()
    generated_thresholds = list((project_dir.parent / "thresholds").glob("thresholds_*.tsv"))
    assert len(generated_thresholds) == 1


def test_generated_project_output_folders_are_discovery_artifacts(tmp_path):
    project_dir = tmp_path / "qupath_dir"
    generated_path = project_dir / "qupath_dir_run_260626_1405" / "tables" / "measurements.csv"

    assert is_qxy_output_artifact(generated_path, project_dir)


def test_latest_timestamped_output_dir_finds_project_run_folders(tmp_path):
    older = tmp_path / "qupath_dir_run_260626_1300"
    newer = tmp_path / "qupath_dir_run_260626_1400"
    check = tmp_path / "qupath_dir_check_260626_1500"
    older.mkdir()
    newer.mkdir()
    check.mkdir()

    assert latest_timestamped_output_dir(tmp_path) == newer
