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
