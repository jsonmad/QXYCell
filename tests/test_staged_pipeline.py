from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def _adata_with_output(tmp_path: Path):
    output_dir = tmp_path / "qxy_outputs_test"
    h5ad_path = output_dir / "h5ad" / "qxycell_test.h5ad"
    obs = pd.DataFrame(
        {"Image": ["image-1"], "Object ID": [1], "CD3_pos": [True]},
        index=["image-1::1"],
    )
    var = pd.DataFrame(
        {"source_measurement_column": ["Cell: CD3 mean"]},
        index=["CD3"],
    )
    adata = ad.AnnData(X=np.array([[2.0]]), obs=obs, var=var)
    adata.uns["qxycell"] = {
        "project_dir": str(tmp_path / "project"),
        "output_dir": str(output_dir),
        "run_dir": str(output_dir),
        "h5ad_path": str(h5ad_path),
        "tables_dir": str(output_dir / "tables"),
    }
    return adata


def _synthetic_project(tmp_path: Path) -> Path:
    project = tmp_path / "project"
    project.mkdir()
    pd.DataFrame(
        {
            "Image": ["image-1", "image-1"],
            "Object ID": [1, 2],
            "Centroid X µm": [5.0, 50.0],
            "Centroid Y µm": [5.0, 50.0],
            "Cell: CD3 mean": [2.0, 8.0],
        }
    ).to_csv(project / "measurements.csv", index=False)
    return project


def _write_annotation(project: Path, label: str) -> None:
    data = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]]],
                },
                "properties": {
                    "objectType": "annotation",
                    "classification": {"name": label},
                },
            }
        ],
    }
    (project / "image-1.geojson").write_text(json.dumps(data), encoding="utf-8")


def _write_classifier(project: Path, threshold: float) -> Path:
    path = project / "CD3.json"
    path.write_text(
        json.dumps(
            {
                "function": {
                    "classifier_fun": "ClassifyByMeasurementFunction",
                    "measurement": "Cell: CD3 mean",
                    "threshold": threshold,
                }
            }
        ),
        encoding="utf-8",
    )
    return path


def _write_threshold_table(project: Path, threshold: float) -> Path:
    path = project / "thresholds.csv"
    pd.DataFrame(
        [
            {
                "marker": "CD3",
                "measurement_column": "Cell: CD3 mean",
                "threshold": threshold,
            }
        ]
    ).to_csv(path, index=False)
    return path


def test_nucleus_classifier_records_compartment_marker_mapping(tmp_path):
    from qxycell.pipeline import import_cells, threshold_from_classifiers

    project = tmp_path / "project"
    project.mkdir()
    pd.DataFrame(
        {
            "Image": ["image-1", "image-1"],
            "Object ID": [1, 2],
            "Centroid X µm": [5.0, 50.0],
            "Centroid Y µm": [5.0, 50.0],
            "Cell: CD3 - Cy5: Median": [2.0, 8.0],
            "Cell: CD3 - Cy5: Mean": [2.0, 8.0],
            "Nucleus: CD3 - Cy5: Median": [2.0, 8.0],
            "Nucleus: CD3 - Cy5: Mean": [2.0, 8.0],
            "Membrane: CD3 - Cy5: Median": [2.0, 8.0],
            "Membrane: CD3 - Cy5: Mean": [2.0, 8.0],
            "Cytoplasm: CD3 - Cy5: Median": [2.0, 8.0],
            "Cytoplasm: CD3 - Cy5: Mean": [2.0, 8.0],
        }
    ).to_csv(project / "measurements.csv", index=False)
    (project / "CD3.json").write_text(
        json.dumps(
            {
                "function": {
                    "classifier_fun": "ClassifyByMeasurementFunction",
                    "measurement": "Nucleus: CD3 - Cy5: Median",
                    "threshold": 5.0,
                }
            }
        ),
        encoding="utf-8",
    )

    adata = import_cells(project, output_dir=tmp_path / "output", verbose=False)
    threshold_from_classifiers(adata, verbose=False)

    assert "CD3-nuc_pos" in adata.obs
    assert "CD3_pos" not in adata.obs
    source_row = adata.var["source_measurement_column"].eq(
        "Nucleus: CD3 - Cy5: Median"
    )
    assert source_row.sum() == 1
    mapped = adata.var.loc[source_row].iloc[0]
    assert mapped["classifier_name"] == "CD3"
    assert mapped["threshold_marker_name"] == "CD3-nuc"
    assert mapped["positivity_column"] == "CD3-nuc_pos"
    assert adata.uns["qxycell_thresholding"]["pos_columns"] == ["CD3-nuc_pos"]


def test_complete_stage_records_columns_and_refreshes_checkpoint(tmp_path):
    from qxycell.stage_state import complete_stage

    adata = _adata_with_output(tmp_path)

    complete_stage(
        adata,
        "thresholds",
        columns=["CD3_pos"],
        details={"source_kind": "classifiers"},
    )

    stage = adata.uns["qxycell"]["stages"]["thresholds"]
    assert stage["status"] == "complete"
    assert stage["columns"] == ["CD3_pos"]
    assert stage["source_kind"] == "classifiers"

    output_dir = Path(adata.uns["qxycell"]["output_dir"])
    saved = ad.read_h5ad(adata.uns["qxycell"]["h5ad_path"])
    assert saved.obs["CD3_pos"].tolist() == [True]
    cells = pd.read_csv(output_dir / "tables" / "cells_obs.csv", index_col=0)
    assert cells["CD3_pos"].tolist() == [True]
    markers = pd.read_csv(output_dir / "tables" / "markers_var.csv", index_col=0)
    assert markers.index.tolist() == ["CD3"]


def test_import_cells_creates_measurement_only_checkpoint(tmp_path):
    from qxycell.pipeline import import_cells

    project = _synthetic_project(tmp_path)
    adata = import_cells(
        project,
        output_dir=tmp_path / "output",
        verbose=False,
    )

    assert adata.shape == (2, 1)
    assert adata.var_names.tolist() == ["CD3 mean"]
    assert adata.obsm["spatial"].tolist() == [[5.0, 5.0], [50.0, 50.0]]
    assert not any(str(column).startswith("annotation__") for column in adata.obs.columns)
    assert adata.uns["qxycell"]["stages"]["measurements"]["status"] == "complete"
    assert Path(adata.uns["qxycell"]["h5ad_path"]).is_file()
    saved = ad.read_h5ad(adata.uns["qxycell"]["h5ad_path"])
    assert saved.uns["qxycell"]["h5ad_path"] == adata.uns["qxycell"]["h5ad_path"]


def test_add_annotations_replaces_previous_geojson_columns(tmp_path):
    from qxycell.pipeline import add_annotations, import_cells

    project = _synthetic_project(tmp_path)
    _write_annotation(project, "Tumor")
    adata = import_cells(project, output_dir=tmp_path / "output", verbose=False)

    add_annotations(adata, pixel_size_um=1.0, verbose=False)
    assert adata.obs["annotation__Tumor"].tolist() == [True, False]

    _write_annotation(project, "Stroma")
    add_annotations(adata, pixel_size_um=1.0, verbose=False)

    assert "annotation__Tumor" not in adata.obs.columns
    assert adata.obs["annotation__Stroma"].tolist() == [True, False]
    saved = ad.read_h5ad(adata.uns["qxycell"]["h5ad_path"])
    assert "annotation__Tumor" not in saved.obs.columns
    assert saved.obs["annotation__Stroma"].tolist() == [True, False]


def test_threshold_from_classifiers_ignores_threshold_tables(tmp_path):
    from qxycell.classifiers import parse_threshold_file
    from qxycell.pipeline import import_cells, threshold_from_classifiers

    project = _synthetic_project(tmp_path)
    _write_classifier(project, 5.0)
    _write_threshold_table(project, 100.0)
    adata = import_cells(project, output_dir=tmp_path / "output", verbose=False)

    summary = threshold_from_classifiers(adata, verbose=False)

    assert adata.obs["CD3_pos"].tolist() == [False, True]
    assert summary["threshold_source_kind"] == "object_classifiers"
    assert summary["threshold_source"].endswith("CD3.json")
    threshold_table = (
        Path(adata.uns["qxycell"]["output_dir"])
        / "thresholds"
        / "classifier_thresholds.tsv"
    )
    assert threshold_table.is_file()
    assert summary["generated_threshold_template"] == str(threshold_table.resolve())
    saved = [item for item in parse_threshold_file(threshold_table) if item.is_simple]
    assert [(item.measurement_column, item.threshold) for item in saved] == [
        ("Cell: CD3 mean", 5.0)
    ]
    assert str(threshold_table.resolve()) in (
        adata.uns["qxycell"]["stages"]["thresholds"]["files"]
    )

    _write_classifier(project, 7.0)
    rerun_summary = threshold_from_classifiers(adata, verbose=False)

    assert rerun_summary["generated_threshold_template"] == str(threshold_table.resolve())
    replaced = [item for item in parse_threshold_file(threshold_table) if item.is_simple]
    assert [(item.measurement_column, item.threshold) for item in replaced] == [
        ("Cell: CD3 mean", 7.0)
    ]
    assert list(threshold_table.parent.glob("classifier_thresholds*.tsv")) == [
        threshold_table
    ]


def test_threshold_from_table_uses_only_named_table_and_replaces_downstream(tmp_path):
    from qxycell.pipeline import (
        import_cells,
        threshold_from_classifiers,
        threshold_from_table,
    )

    project = _synthetic_project(tmp_path)
    _write_classifier(project, 5.0)
    table = _write_threshold_table(project, 1.0)
    adata = import_cells(project, output_dir=tmp_path / "output", verbose=False)
    threshold_from_classifiers(adata, verbose=False)
    adata.obs["celltype"] = ["Low", "High"]
    adata.uns["qxycell"]["stages"]["celltypes"] = {
        "status": "complete",
        "columns": ["celltype"],
        "files": [],
    }
    adata.uns["qxycell"]["celltyping_applied"] = True
    adata.uns["qxycell_celltyping"] = {"celltype_column": "celltype"}

    summary = threshold_from_table(adata, table, verbose=False)

    assert adata.obs["CD3_pos"].tolist() == [True, True]
    assert "celltype" not in adata.obs.columns
    assert adata.uns["qxycell"]["stages"]["celltypes"]["status"] == "stale"
    assert adata.uns["qxycell"]["celltyping_applied"] is False
    assert "qxycell_celltyping" not in adata.uns
    assert summary["threshold_source_kind"] == "manual_threshold_file"
    assert summary["threshold_source"] == str(table.resolve())
    cells = pd.read_csv(
        Path(adata.uns["qxycell"]["tables_dir"]) / "cells_obs.csv",
        index_col=0,
    )
    assert cells["CD3_pos"].tolist() == [True, True]
    assert "celltype" not in cells.columns


def test_celltype_prompt_replaces_current_prompt_and_preserves_yaml(tmp_path):
    from qxycell.pipeline import import_cells, threshold_from_classifiers
    from qxycell.prompts import celltype_prompt

    project = _synthetic_project(tmp_path)
    _write_classifier(project, 5.0)
    adata = import_cells(project, output_dir=tmp_path / "output", verbose=False)
    threshold_from_classifiers(adata, verbose=False)
    celltype_dir = Path(adata.uns["qxycell"]["output_dir"]) / "celltype"
    celltype_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = celltype_dir / "celltype_logic.yaml"
    yaml_path.write_text("rules: []\nfeatures: {}\nderived_features: {}\n", encoding="utf-8")

    celltype_prompt(adata, context="first context", print_prompt=False)
    prompt_path = celltype_dir / "current_prompt.txt"
    assert prompt_path.is_file()
    assert "first context" in prompt_path.read_text(encoding="utf-8")

    celltype_prompt(adata, context="second context", print_prompt=False)

    assert "second context" in prompt_path.read_text(encoding="utf-8")
    assert "first context" not in prompt_path.read_text(encoding="utf-8")
    assert yaml_path.read_text(encoding="utf-8").startswith("rules: []")
    prompt_stage = adata.uns["qxycell"]["stages"]["celltype_prompt"]
    assert prompt_stage["status"] == "complete"
    assert prompt_stage["files"] == [str(prompt_path.resolve())]


def test_celltype_rerun_removes_obsolete_feature_columns(tmp_path):
    from qxycell.celltyping import apply_celltypes
    from qxycell.pipeline import import_cells, threshold_from_classifiers

    project = _synthetic_project(tmp_path)
    _write_classifier(project, 5.0)
    adata = import_cells(project, output_dir=tmp_path / "output", verbose=False)
    threshold_from_classifiers(adata, verbose=False)
    first_logic = {
        "rules": [{"name": "T_cell", "positive": ["CD3"]}],
        "features": {"Old_feature": {"positive": ["CD3"]}},
        "derived_features": {},
    }
    second_logic = {
        "rules": [{"name": "T_cell", "positive": ["CD3"]}],
        "features": {"New_feature": {"positive": ["CD3"]}},
        "derived_features": {},
    }

    apply_celltypes(adata, first_logic, verbose=False)
    assert "Old_feature" in adata.obs.columns
    apply_celltypes(adata, second_logic, verbose=False)

    assert "Old_feature" not in adata.obs.columns
    assert adata.obs["New_feature"].tolist() == [0, 1]
    counts_path = Path(adata.uns["qxycell"]["tables_dir"]) / "celltype_counts.csv"
    assert counts_path.is_file()
    saved = ad.read_h5ad(adata.uns["qxycell"]["h5ad_path"])
    assert "Old_feature" not in saved.obs.columns
    assert saved.obs["New_feature"].tolist() == [0, 1]


def test_public_api_exposes_staged_functions():
    import qxycell as qxy

    assert callable(qxy.import_cells)
    assert not hasattr(qxy, "import_measurements")
    assert not hasattr(qxy, "run")
    assert callable(qxy.add_annotations)
    assert callable(qxy.threshold_from_classifiers)
    assert callable(qxy.threshold_from_table)
