import json

import pandas as pd
import pytest

from qxycell.markers import marker_name_from_classifier


@pytest.mark.parametrize(
    ("measurement", "expected"),
    [
        ("Cell: CD3 - Cy5: Median", "CD3"),
        ("Nucleus: CD3 - Cy5: Median", "CD3-nuc"),
        ("Membrane: CD3 - Cy5: Mean", "CD3-mem"),
        ("Cytoplasm: CD3 - Cy5: Mean", "CD3-cyto"),
        ("CD3 intensity", "CD3"),
        (None, "CD3"),
    ],
)
def test_marker_name_from_classifier_adds_compartment_suffix(measurement, expected):
    assert marker_name_from_classifier("CD3", measurement) == expected


def test_marker_name_from_classifier_is_idempotent():
    assert (
        marker_name_from_classifier("CD3-nuc", "Nucleus: CD3 - Cy5: Median")
        == "CD3-nuc"
    )


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
