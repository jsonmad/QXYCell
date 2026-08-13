"""Marker naming helpers."""

from __future__ import annotations

import re


def marker_name_from_classifier_name(name: str) -> str:
    """Normalize a classifier filename stem into a valid marker variable name."""

    cleaned = re.sub(r"[^0-9A-Za-z]+", "_", str(name).strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned or "marker"

