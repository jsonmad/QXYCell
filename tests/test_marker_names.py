import json
from pathlib import Path

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


def test_marker_allocators_avoid_existing_numeric_suffix_collisions():
    from qxycell.pipeline import (
        _unique_marker_names,
        _unique_marker_names_for_measurement_columns,
    )
    from qxycell.types import ClassifierDefinition

    measurement_columns = [
        "Cell: first mean",
        "Cell: second mean",
        "Cell: third mean",
    ]
    classifiers = [
        ClassifierDefinition(
            path=Path(f"classifier-{index}.json"),
            name=name,
            measurement_column=measurement,
            threshold=1.0,
            is_simple=True,
        )
        for index, (name, measurement) in enumerate(
            zip(
                ["CD3", "CD3", "CD3_2"],
                measurement_columns,
                strict=True,
            )
        )
    ]
    expected = {0: "CD3", 1: "CD3_2", 2: "CD3_2_2"}

    assert _unique_marker_names([[item] for item in classifiers]) == expected
    assert (
        _unique_marker_names_for_measurement_columns(
            measurement_columns,
            classifiers,
        )
        == expected
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


def test_classifier_table_round_trip_preserves_exact_measurement_ownership(tmp_path):
    from qxycell.pipeline import (
        import_cells,
        threshold_from_classifiers,
        threshold_from_table,
    )

    mean_measurement = "Nucleus: CD3 - Cy5: Mean"
    median_measurement = "Nucleus: CD3 - Cy5: Median"
    project = tmp_path / "project"
    project.mkdir()
    pd.DataFrame(
        {
            "Image": ["image-1", "image-1"],
            "Object ID": [1, 2],
            "Centroid X µm": [5.0, 50.0],
            "Centroid Y µm": [5.0, 50.0],
            mean_measurement: [1.0, 9.0],
            median_measurement: [9.0, 1.0],
        }
    ).to_csv(project / "measurements.csv", index=False)
    for filename, measurement, threshold in (
        ("CD3!.json", median_measurement, 7.0),
        ("CD3.json", mean_measurement, 3.0),
    ):
        (project / filename).write_text(
            json.dumps(
                {
                    "function": {
                        "classifier_fun": "ClassifyByMeasurementFunction",
                        "measurement": measurement,
                        "threshold": threshold,
                    }
                }
            ),
            encoding="utf-8",
        )

    json_adata = import_cells(
        project,
        output_dir=tmp_path / "json-output",
        verbose=False,
    )
    json_summary = threshold_from_classifiers(json_adata, verbose=False)
    json_ownership = {
        marker: source
        for marker, source in zip(
            json_adata.var["threshold_marker_name"],
            json_adata.var["source_measurement_column"],
            strict=True,
        )
        if marker
    }

    table_adata = import_cells(
        project,
        output_dir=tmp_path / "table-output",
        verbose=False,
    )
    threshold_from_table(
        table_adata,
        json_summary["generated_threshold_template"],
        verbose=False,
    )
    table_ownership = {
        marker: source
        for marker, source in zip(
            table_adata.var["threshold_marker_name"],
            table_adata.var["source_measurement_column"],
            strict=True,
        )
        if marker
    }

    assert json_ownership["CD3-nuc"] == mean_measurement
    assert table_ownership == json_ownership
