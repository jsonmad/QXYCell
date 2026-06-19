import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from qxycell.checks import check
from qxycell.geojson import load_cell_polygons, summarize_geojson_file
from qxycell.pipeline import _apply_annotations, _apply_cell_polygons, run
from qxycell.tma import assign_tma_cores


def _feature(*, object_type, name, coords, classification=None):
    properties = {}
    if object_type is not None:
        properties["objectType"] = object_type
    if name is not None:
        properties["name"] = name
    if classification is not None:
        properties["classification"] = {"name": classification}
    return {
        "type": "Feature",
        "properties": properties,
        "geometry": {"type": "Polygon", "coordinates": [coords]},
    }


def _write_geojson(path, features):
    path.write_text(
        json.dumps({"type": "FeatureCollection", "features": features}),
        encoding="utf-8",
    )
    return path


def _adata_for_img(coords):
    obs = pd.DataFrame(
        {
            "Image": ["img.ome.tiff"] * len(coords),
            "Object ID": [f"cell_{i}" for i in range(len(coords))],
            "Centroid X µm": [x for x, _ in coords],
            "Centroid Y µm": [y for _, y in coords],
        },
        index=[f"img.ome.tiff::cell_{i}" for i in range(len(coords))],
    )
    adata = ad.AnnData(X=np.zeros((len(coords), 1)), obs=obs)
    adata.obsm["spatial"] = obs[["Centroid X µm", "Centroid Y µm"]].to_numpy()
    return adata


def _write_minimal_run_project(project_dir):
    project_dir.mkdir()
    pd.DataFrame(
        {
            "Image": ["img.ome.tiff", "img.ome.tiff"],
            "Object ID": ["cell_0", "cell_1"],
            "Centroid X µm": [5, 25],
            "Centroid Y µm": [5, 25],
            "Marker: Mean": [10.0, 1.0],
        }
    ).to_csv(project_dir / "detections.tsv", sep="\t", index=False)
    (project_dir / "marker.json").write_text(
        json.dumps(
            {
                "function": {
                    "classifier_fun": "ClassifyByMeasurementFunction",
                    "measurement": "Marker: Mean",
                    "threshold": 5,
                }
            }
        ),
        encoding="utf-8",
    )


def _write_minimal_manual_threshold_project(project_dir):
    project_dir.mkdir()
    pd.DataFrame(
        {
            "Image": ["img_a.ome.tiff", "img_a.ome.tiff", "img_b.ome.tiff", "img_b.ome.tiff"],
            "Object ID": ["cell_0", "cell_1", "cell_2", "cell_3"],
            "Centroid X µm": [5, 25, 5, 25],
            "Centroid Y µm": [5, 25, 5, 25],
            "Marker: Mean": [10.0, 1.0, 10.0, 1.0],
        }
    ).to_csv(project_dir / "detections.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "marker": ["Marker"],
            "measurement_column": ["Marker: Mean"],
            "img_a.ome.tiff": [5],
            "img_b.ome.tiff": [11],
        }
    ).to_csv(project_dir / "thresholds_260615-1234.tsv", sep="\t", index=False)


def _write_template_project(project_dir):
    project_dir.mkdir()
    pd.DataFrame(
        {
            "Image": ["img_a.ome.tiff", "img_b.ome.tiff"],
            "Object ID": ["cell_0", "cell_1"],
            "Centroid X µm": [5, 25],
            "Centroid Y µm": [5, 25],
            "Cell: CD8 - Cy5: Mean": [8.0, 2.0],
            "Nucleus: DAPI: Mean": [100.0, 50.0],
            "Marker: Mean": [10.0, 1.0],
            "Marker: Median": [9.0, 1.5],
            "Marker: Area": [100.0, 200.0],
        }
    ).to_csv(project_dir / "detections.tsv", sep="\t", index=False)
    (project_dir / "marker.json").write_text(
        json.dumps(
            {
                "function": {
                    "classifier_fun": "ClassifyByMeasurementFunction",
                    "measurement": "Marker: Mean",
                    "threshold": 5,
                }
            }
        ),
        encoding="utf-8",
    )


def test_run_annotation_mapper_ignores_tma_core_features(tmp_path):
    geojson_path = _write_geojson(
        tmp_path / "img.geojson",
        [
            _feature(
                object_type="annotation",
                name="Region*",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="tmaCore",
                name="core_1",
                coords=[[0, 0], [20, 0], [20, 20], [0, 20], [0, 0]],
            ),
        ],
    )
    adata = _adata_for_img([(5, 5)])

    conflicts = _apply_annotations(
        adata,
        [summarize_geojson_file(geojson_path)],
        pixel_size_um=1.0,
    )

    assert conflicts == []
    assert bool(adata.obs.loc["img.ome.tiff::cell_0", "annotation__Region"])
    assert "tma_core" not in adata.obs.columns


def test_run_detects_tma_cores_without_auto_assigning(tmp_path):
    project_dir = tmp_path / "project"
    _write_minimal_run_project(project_dir)
    _write_geojson(
        project_dir / "img.geojson",
        [
            _feature(
                object_type="tmaCore",
                name="core_1",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            )
        ],
    )

    adata = run(
        project_dir,
        output_dir=tmp_path / "out",
        pixel_size_um=1.0,
        verbose=False,
    )

    assert "Xµm" in adata.obs.columns
    assert "Yµm" in adata.obs.columns
    assert "Centroid X µm" not in adata.obs.columns
    assert "Centroid Y µm" not in adata.obs.columns
    assert adata.obs["Xµm"].tolist() == [5, 25]
    assert adata.obs["Yµm"].tolist() == [5, 25]
    np.testing.assert_array_equal(adata.obsm["spatial"], np.array([[5, 5], [25, 25]]))
    assert "tma_core" not in adata.obs.columns
    assert adata.uns["qxycell"]["n_tma_core_features"] == 1
    assert adata.uns["qxycell"]["tma_assignment"] is None
    assert "qxycell_tma" not in adata.uns


def test_check_writes_manual_threshold_template(tmp_path):
    project_dir = tmp_path / "project"
    _write_template_project(project_dir)

    report = check(project_dir, output_dir=tmp_path / "qxy_outputs_260615-1234")

    template_path = report.output_dir / "tables" / "thresholds_260615-1234.tsv"
    assert template_path.exists()
    template = pd.read_csv(template_path, sep="\t")
    assert template.columns.tolist() == [
        "compartment",
        "marker",
        "measurement_column",
        "img_a.ome.tiff",
        "img_b.ome.tiff",
    ]
    assert template["compartment"].fillna("").tolist() == ["Cell", "", "", "Nucleus"]
    assert template["marker"].tolist() == ["CD8", "Marker", "Marker", "DAPI"]
    assert template["measurement_column"].tolist() == [
        "Cell: CD8 - Cy5: Mean",
        "Marker: Mean",
        "Marker: Median",
        "Nucleus: DAPI: Mean",
    ]
    assert template["img_a.ome.tiff"].tolist()[1] == 5.0
    assert template["img_b.ome.tiff"].tolist()[1] == 5.0
    assert template[["img_a.ome.tiff", "img_b.ome.tiff"]].iloc[2].isna().tolist() == [
        True,
        True,
    ]


def test_run_uses_manual_threshold_tsv_when_simple_json_missing(tmp_path):
    project_dir = tmp_path / "project"
    _write_minimal_manual_threshold_project(project_dir)

    adata = run(
        project_dir,
        output_dir=tmp_path / "out",
        pixel_size_um=1.0,
        verbose=False,
    )

    assert adata.var_names.tolist() == ["Marker"]
    assert adata.var.loc["Marker", "source_measurement_column"] == "Marker: Mean"
    assert adata.var.loc["Marker", "threshold"] == "per_image"
    assert adata.obs["Marker_pos"].tolist() == [1, 0, 0, 0]
    assert adata.uns["qxycell"]["n_simple_classifiers"] == 2


def test_run_uses_newest_timestamped_threshold_file(tmp_path):
    project_dir = tmp_path / "project"
    _write_minimal_manual_threshold_project(project_dir)
    older = project_dir / "thresholds_260615-1234.tsv"
    newer = project_dir / "thresholds_260615-1300.tsv"
    pd.DataFrame(
        {
            "marker": ["Marker"],
            "measurement_column": ["Marker: Mean"],
            "img_a.ome.tiff": [11],
            "img_b.ome.tiff": [11],
        }
    ).to_csv(newer, sep="\t", index=False)
    os.utime(older, (1_700_000_000, 1_700_000_000))
    os.utime(newer, (1_700_000_100, 1_700_000_100))

    report = check(project_dir, output_dir=tmp_path / "out")
    ignored = [
        message
        for message in report.messages
        if message.code == "classifiers.threshold_file_ignored"
    ]
    assert len(ignored) == 1
    assert ignored[0].path == str(older.resolve())

    adata = run(
        project_dir,
        output_dir=tmp_path / "run_out",
        pixel_size_um=1.0,
        verbose=False,
    )

    assert adata.obs["Marker_pos"].tolist() == [0, 0, 0, 0]


def test_run_sample_annotations_collapse_to_one_sample_column(tmp_path):
    project_dir = tmp_path / "project"
    _write_minimal_run_project(project_dir)
    _write_geojson(
        project_dir / "img.geojson",
        [
            _feature(
                object_type="annotation",
                name="Region*",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="annotation",
                name="Ignore*",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="annotation",
                name="Sample-A",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="annotation",
                name="Sample-B",
                coords=[[0, 0], [30, 0], [30, 30], [0, 30], [0, 0]],
            ),
        ],
    )

    with pytest.warns(UserWarning, match="more than one sample annotation"):
        adata = run(
            project_dir,
            output_dir=tmp_path / "out",
            pixel_size_um=1.0,
            verbose=False,
        )

    assert bool(adata.obs.loc["img.ome.tiff::cell_0", "annotation__Region"])
    assert bool(adata.obs.loc["img.ome.tiff::cell_0", "annotation__Ignore"])
    assert "annotation__Sample_A" not in adata.obs.columns
    assert "sample_annotation__Sample_A" not in adata.obs.columns
    assert adata.obs["Sample"].astype(str).tolist() == ["Ambiguous", "Sample-B"]
    assert adata.uns["qxycell_sample_annotations"]["n_conflicting_cells"] == 1
    assert adata.uns["qxycell"]["n_annotation_conflicts"] == 0


def test_run_annotation_mapper_uses_classification_and_name_labels(tmp_path):
    geojson_path = _write_geojson(
        tmp_path / "img.geojson",
        [
            _feature(
                object_type="annotation",
                name="artifacts",
                classification="Ignore*",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
        ],
    )
    adata = _adata_for_img([(5, 5), (20, 20)])

    _apply_annotations(
        adata,
        [summarize_geojson_file(geojson_path)],
        pixel_size_um=1.0,
    )

    assert adata.obs["annotation__Ignore"].tolist() == [True, False]
    assert adata.obs["annotation__artifacts"].tolist() == [True, False]
    assert adata.uns["qxycell_annotation_labels"]["annotation__Ignore"] == "Ignore*"
    assert adata.uns["qxycell_annotation_labels"]["annotation__artifacts"] == "artifacts"


def test_run_cell_polygons_match_qupath_cells_suffix_by_object_id(tmp_path):
    geojson_path = _write_geojson(
        tmp_path / "img-cells.geojson",
        [
            {
                **_feature(
                    object_type="cell",
                    name=None,
                    coords=[[0, 0], [2, 0], [2, 2], [0, 2], [0, 0]],
                ),
                "id": "cell_0",
            },
            {
                **_feature(
                    object_type="cell",
                    name=None,
                    coords=[[10, 10], [12, 10], [12, 12], [10, 12], [10, 10]],
                ),
                "id": "cell_1",
            },
        ],
    )
    adata = _adata_for_img([(1, 1), (11, 11), (30, 30)])

    n_matched = _apply_cell_polygons(
        adata,
        [summarize_geojson_file(geojson_path)],
        pixel_size_um=1.0,
    )

    assert n_matched == 2
    assert adata.obs["cell_polygon_wkt"].str.startswith("POLYGON").tolist() == [
        True,
        True,
        False,
    ]
    assert "cell_polygons" not in adata.obsm

    output_path = tmp_path / "cells.h5ad"
    adata.write_h5ad(output_path)
    reloaded = ad.read_h5ad(output_path)
    assert reloaded.obs["cell_polygon_wkt"].str.startswith("POLYGON").tolist() == [
        True,
        True,
        False,
    ]


def test_load_cell_polygons_defaults_to_obs_wkt_column(tmp_path):
    _write_geojson(
        tmp_path / "img-cells.geojson",
        [
            {
                **_feature(
                    object_type="cell",
                    name=None,
                    coords=[[0, 0], [2, 0], [2, 2], [0, 2], [0, 0]],
                ),
                "id": "cell_0",
            }
        ],
    )
    adata = _adata_for_img([(1, 1), (11, 11)])

    n_matched = load_cell_polygons(
        adata,
        tmp_path,
        pixel_size_um=1.0,
        verbose=False,
    )

    assert n_matched == 1
    assert adata.obs["cell_polygon_wkt"].str.startswith("POLYGON").tolist() == [
        True,
        False,
    ]
    assert "cell_polygons" not in adata.obsm


def test_assign_tma_cores_leaves_overlap_cells_unassigned(tmp_path):
    geojson_path = _write_geojson(
        tmp_path / "img.geojson",
        [
            _feature(
                object_type="tmaCore",
                name="core_1",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="tmaCore",
                name="core_2",
                coords=[[5, 0], [15, 0], [15, 10], [5, 10], [5, 0]],
            ),
        ],
    )
    adata = _adata_for_img([(2, 5), (7, 5), (12, 5)])

    summary = assign_tma_cores(
        adata,
        geojson_path,
        pixel_size_um=1.0,
        output_dir=tmp_path / "out",
        verbose=False,
    )

    assert adata.obs["tma_core"].tolist() == ["core_1", "Unassigned", "core_2"]
    assert summary["n_overlapping_core_pairs"] == 1
    assert summary["n_ambiguous_overlap_cells"] == 1
    assert summary["n_assigned_cells"] == 2


def test_assign_tma_cores_ignores_non_tma_annotation_features_by_default(tmp_path):
    geojson_path = _write_geojson(
        tmp_path / "img.geojson",
        [
            _feature(
                object_type="annotation",
                name="Ignore*",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="tmaCore",
                name="core_1",
                coords=[[20, 20], [30, 20], [30, 30], [20, 30], [20, 20]],
            ),
        ],
    )
    adata = _adata_for_img([(5, 5), (25, 25)])

    summary = assign_tma_cores(
        adata,
        geojson_path,
        pixel_size_um=1.0,
        output_dir=tmp_path / "out",
        verbose=False,
    )

    assert summary["n_cores"] == 1
    assert adata.obs["tma_core"].tolist() == ["Unassigned", "core_1"]


def test_assign_tma_cores_includes_boundary_points(tmp_path):
    geojson_path = _write_geojson(
        tmp_path / "img.geojson",
        [
            _feature(
                object_type="tmaCore",
                name="core_1",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            )
        ],
    )
    adata = _adata_for_img([(0, 5)])

    assign_tma_cores(
        adata,
        geojson_path,
        pixel_size_um=1.0,
        output_dir=tmp_path / "out",
        verbose=False,
    )

    assert adata.obs["tma_core"].tolist() == ["core_1"]


def test_check_report_lists_tma_core_labels_under_tma(tmp_path):
    project_dir = tmp_path / "project"
    _write_minimal_run_project(project_dir)
    _write_geojson(
        project_dir / "img.geojson",
        [
            _feature(
                object_type="tmaCore",
                name="A-1",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="annotation",
                name="sample1",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="annotation",
                name="Ignore*",
                coords=[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
            ),
            _feature(
                object_type="cell",
                name="Ignore*",
                coords=[[0, 0], [2, 0], [2, 2], [0, 2], [0, 0]],
            ),
            _feature(
                object_type="cell",
                name="Ignore*",
                coords=[[3, 3], [4, 3], [4, 4], [3, 4], [3, 3]],
            ),
        ],
    )

    report = check(project_dir, output_dir=tmp_path / "check")
    report_text = report.report_path.read_text(encoding="utf-8")

    assert "Sample   : sample1 (1)" in report_text
    assert "Ignore   : Ignore* (1)" in report_text
    assert "TMA      : A-1 (1)" in report_text
    assert "Other    : none" in report_text
    assert "Cell labels: Ignore* (2)" in report_text
