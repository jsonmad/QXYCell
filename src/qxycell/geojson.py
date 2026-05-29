"""GeoJSON discovery and summaries."""

from __future__ import annotations

import collections
import json
from pathlib import Path
from typing import Any

from qxycell.types import GeoJsonFile, Message


def discover_geojson_files(project_dir: str | Path) -> list[Path]:
    """Find exported QuPath GeoJSON files."""

    root = Path(project_dir).expanduser().resolve()
    files = [path for path in root.rglob("*.geojson") if not path.name.startswith(".")]
    return sorted(dict.fromkeys(files))


def _classification_name(properties: dict[str, Any]) -> str:
    classification = properties.get("classification")
    if isinstance(classification, dict):
        label = classification.get("name")
        if label:
            return str(label)
    elif classification is not None:
        return str(classification)

    name = properties.get("name")
    if name:
        return str(name)
    return "Unclassified"


def summarize_geojson_file(path: str | Path) -> GeoJsonFile:
    """Summarize feature counts, object types, class names, and object names."""

    path = Path(path).expanduser().resolve()
    try:
        data = json.loads(path.read_text(errors="replace"))
        features = data.get("features", []) if isinstance(data, dict) else []
        object_type_counts: collections.Counter[str] = collections.Counter()
        class_counts: collections.Counter[str] = collections.Counter()
        name_counts: collections.Counter[str] = collections.Counter()
        for feature in features:
            properties = feature.get("properties", {}) if isinstance(feature, dict) else {}
            if not isinstance(properties, dict):
                properties = {}
            object_type_counts[str(properties.get("objectType") or "unknown")] += 1
            class_counts[_classification_name(properties)] += 1
            name_counts[str(properties.get("name") or "")] += 1
        return GeoJsonFile(
            path=path,
            n_features=len(features),
            object_type_counts=dict(object_type_counts),
            class_counts=dict(class_counts),
            name_counts=dict(name_counts),
        )
    except Exception as exc:
        return GeoJsonFile(path=path, n_features=None, readable=False, error=str(exc))


def summarize_geojson_files(paths: list[Path]) -> list[GeoJsonFile]:
    """Summarize all GeoJSON files."""

    return [summarize_geojson_file(path) for path in paths]


def validate_geojson_files(files: list[GeoJsonFile]) -> list[Message]:
    """Validate GeoJSON readability."""

    messages: list[Message] = []
    if not files:
        messages.append(
            Message(
                level="warning",
                code="geojson.missing",
                message="No GeoJSON annotation files were found.",
            )
        )
    for file in files:
        if not file.readable:
            messages.append(
                Message(
                    level="error",
                    code="geojson.unreadable",
                    message=f"GeoJSON file could not be read: {file.error}",
                    path=str(file.path),
                )
            )
    return messages
