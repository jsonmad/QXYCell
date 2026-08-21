"""Marker naming helpers."""

from __future__ import annotations

import re


COMPARTMENT_SUFFIXES = {
    "Nucleus": "nuc",
    "Membrane": "mem",
    "Cytoplasm": "cyto",
}


def marker_name_from_classifier_name(name: str) -> str:
    """Normalize a classifier filename stem into a valid marker variable name."""

    cleaned = re.sub(r"[^0-9A-Za-z]+", "_", str(name).strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned or "marker"


def marker_name_from_classifier(name: str, measurement_column: str | None) -> str:
    """Return the canonical marker name for a classifier measurement."""

    from qxycell.classifiers import compartment_from_measurement_column

    compartment = compartment_from_measurement_column(measurement_column or "")
    suffix = COMPARTMENT_SUFFIXES.get(compartment)
    raw_name = str(name).strip()
    if suffix is None:
        return marker_name_from_classifier_name(raw_name)
    raw_name = re.sub(rf"[-_]{re.escape(suffix)}$", "", raw_name, flags=re.IGNORECASE)
    marker = marker_name_from_classifier_name(raw_name)
    return f"{marker}-{suffix}"
