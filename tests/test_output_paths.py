import json
import re

import anndata as ad
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
    assert (output_dir / "h5ad" / "qxycell.h5ad").exists()
    assert not (output_dir / "run").exists()
    assert adata.uns["qxycell"]["check_output_dir"] is None
    assert adata.uns["qxycell"]["check_report_txt"] is None
    assert adata.uns["qxycell"]["validation_ok"] is True
    assert "check_ok" not in adata.uns["qxycell"]
    assert not any(project_dir.parent.glob("qupath_dir_check_*"))
    assert not (output_dir / "check_report.txt").exists()
    generated_thresholds = list((output_dir / "thresholds").glob("thresholds_*.tsv"))
    assert len(generated_thresholds) == 1


def test_workflow_accepts_and_applies_explicit_threshold_file(tmp_path):
    project_dir = tmp_path / "qupath_dir"
    output_dir = tmp_path / "outputs"
    _write_minimal_project(project_dir)
    threshold_file = tmp_path / "reviewed_thresholds.tsv"
    pd.DataFrame(
        {
            "marker": ["CD3"],
            "measurement_column": ["Cell: CD3: Mean"],
            "sample_A.tif": [0.9],
        }
    ).to_csv(threshold_file, sep="\t", index=False)

    adata = qxy.workflow(
        project_dir,
        output_dir=output_dir,
        threshold_file=threshold_file,
        apply_thresholds=True,
        remove_ignore_cells=False,
        make_qc=False,
        make_plots=False,
        verbose=False,
    )

    assert adata.obs["CD3_pos"].tolist() == [False]
    assert adata.uns["qxycell_workflow"]["threshold_file"] == str(threshold_file)
    assert adata.uns["qxycell_workflow"]["thresholds_applied"] is True

    reloaded = ad.read_h5ad(output_dir / "h5ad" / "qxycell.h5ad")
    assert reloaded.uns["qxycell_workflow"] == adata.uns["qxycell_workflow"]


def test_generated_project_output_folders_are_discovery_artifacts(tmp_path):
    project_dir = tmp_path / "qupath_dir"
    generated_path = project_dir / "qupath_dir_run_260626_1405" / "tables" / "measurements.csv"

    assert is_qxy_output_artifact(generated_path, project_dir)


def test_explicit_timestamped_output_folders_are_discovery_artifacts(tmp_path):
    project_dir = tmp_path / "qupath_dir"
    generated_path = (
        project_dir
        / "outputs_260713-141927"
        / "thresholds"
        / "thresholds_260713-1419.tsv"
    )

    assert is_qxy_output_artifact(generated_path, project_dir)


def test_latest_timestamped_output_dir_finds_project_run_folders(tmp_path):
    older = tmp_path / "qupath_dir_run_260626_1300"
    newer = tmp_path / "qupath_dir_run_260626_1400"
    check = tmp_path / "qupath_dir_check_260626_1500"
    older.mkdir()
    newer.mkdir()
    check.mkdir()

    assert latest_timestamped_output_dir(tmp_path) == newer
