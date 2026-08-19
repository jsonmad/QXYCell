"""GeoJSON discovery, summaries, and cell polygon loading."""

from __future__ import annotations

import collections
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

from qxycell.discovery import is_qxy_output_artifact
from qxycell.types import GeoJsonFile, Message

if TYPE_CHECKING:
    pass


def discover_geojson_files(project_dir: str | Path) -> list[Path]:
    """Find exported QuPath GeoJSON files."""

    root = Path(project_dir).expanduser().resolve()
    files = [
        path
        for path in root.rglob("*.geojson")
        if not path.name.startswith(".") and not is_qxy_output_artifact(path, root)
    ]
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


def _label_values(properties: dict[str, Any]) -> list[str]:
    labels = [_classification_name(properties)]
    name = properties.get("name")
    if name:
        label = str(name)
        if label not in labels:
            labels.append(label)
    return labels


def summarize_geojson_file(path: str | Path) -> GeoJsonFile:
    """Summarize feature counts, object types, class names, and object names."""

    path = Path(path).expanduser().resolve()
    try:
        data = json.loads(path.read_text(errors="replace"))
        features = data.get("features", []) if isinstance(data, dict) else []
        object_type_counts: collections.Counter[str] = collections.Counter()
        class_counts: collections.Counter[str] = collections.Counter()
        name_counts: collections.Counter[str] = collections.Counter()
        labels_by_object_type: dict[str, collections.Counter[str]] = {}
        n_supported_features = 0
        for feature in features:
            properties = feature.get("properties", {}) if isinstance(feature, dict) else {}
            if not isinstance(properties, dict):
                properties = {}
            object_type = str(properties.get("objectType") or "unknown")
            if object_type.lower() not in {"annotation", "cell"}:
                continue
            n_supported_features += 1
            object_type_counts[object_type] += 1
            class_counts[_classification_name(properties)] += 1
            name_counts[str(properties.get("name") or "")] += 1
            labels_by_object_type.setdefault(object_type, collections.Counter())
            for label in _label_values(properties):
                labels_by_object_type[object_type][label] += 1
        return GeoJsonFile(
            path=path,
            n_features=n_supported_features,
            object_type_counts=dict(object_type_counts),
            class_counts=dict(class_counts),
            name_counts=dict(name_counts),
            labels_by_object_type={
                object_type: dict(labels)
                for object_type, labels in labels_by_object_type.items()
            },
        )
    except Exception as exc:
        return GeoJsonFile(path=path, n_features=None, readable=False, error=str(exc))


def summarize_geojson_files(paths: list[Path]) -> list[GeoJsonFile]:
    """Summarize all GeoJSON files."""

    return [summarize_geojson_file(path) for path in paths]


def load_cell_polygons(
    adata: "Any",
    project_dir: str | Path,
    *,
    object_id_col: str = "Object ID",
    pixel_size_um: float = 0.28,
    obs_key: str = "cell_polygon_wkt",
    obsm_key: str | None = None,
    verbose: bool = True,
) -> int:
    """Load cell segmentation polygons from all GeoJSON files in a QuPath project folder.

    Scans every ``.geojson`` file in ``project_dir``, reads all features that have an
    ``Object ID`` property and a valid polygon geometry, and matches them to rows in
    ``adata.obs`` by Object ID alone (no filename-to-image matching required).

    Matched geometries are converted to micron coordinates (scaled by ``pixel_size_um``)
    and stored as WKT strings in ``adata.obs[obs_key]``. Unmatched cells receive an
    empty string.

    Parameters
    ----------
    adata:
        AnnData object produced by ``qxy.import_cells()``.
    project_dir:
        Path to the QuPath project folder containing ``.geojson`` files.
    object_id_col:
        Column in ``adata.obs`` containing QuPath Object IDs (default ``"Object ID"``).
    pixel_size_um:
        Pixel size in microns used to scale polygon coordinates to microns
        (default ``0.28``). Must match the value used in ``qxy.add_annotations()``.
    obs_key:
        Column under which WKT strings are stored in ``adata.obs``
        (default ``"cell_polygon_wkt"``).
    obsm_key:
        Optional legacy key under which WKT strings are also stored in
        ``adata.obsm``. Leave as ``None`` for H5AD-friendly string storage in
        ``adata.obs`` only.
    verbose:
        Print a summary on completion (default ``True``).

    Returns
    -------
    int
        Number of cells successfully matched to a polygon.

    Examples
    --------
    >>> n = qxy.load_cell_polygons(adata, "/path/to/qupath_project")
    >>> from shapely import wkt
    >>> geom = wkt.loads(adata.obs["cell_polygon_wkt"].iloc[0])

    Convert to a GeoDataFrame for spatial analysis::

        import geopandas as gpd
        gdf = gpd.GeoDataFrame(
            adata.obs,
            geometry=gpd.GeoSeries.from_wkt(adata.obs["cell_polygon_wkt"]),
        )
    """
    try:
        import numpy as _np
        from shapely import affinity as _affinity
        from shapely.geometry import shape as _shape
    except ImportError as exc:
        raise ImportError(
            "load_cell_polygons requires shapely. "
            "Install it with: pip install shapely"
        ) from exc

    project_path = Path(project_dir).expanduser().resolve()
    geojson_paths = discover_geojson_files(project_path)

    # Build lookup: object_id → shapely geometry
    # Object IDs in QuPath are UUIDs and are globally unique across images,
    # so we match on Object ID alone without needing to know which image each
    # GeoJSON belongs to.
    geom_by_id: dict[str, Any] = {}
    n_features_read = 0
    n_skipped = 0

    for path in geojson_paths:
        try:
            data = json.loads(path.read_text(errors="replace"))
        except Exception:
            continue
        features = data.get("features", []) if isinstance(data, dict) else []

        for feature in features:
            if not isinstance(feature, dict):
                continue
            properties = feature.get("properties") or {}
            if not isinstance(properties, dict):
                properties = {}
            if str(properties.get("objectType") or "").lower() != "cell":
                continue

            # QuPath stores the Object ID as the top-level GeoJSON feature "id"
            # field, not inside "properties". Fall back to property keys for
            # other exporters.
            obj_id = (
                feature.get("id")
                or properties.get("objectID")
                or properties.get("id")
                or properties.get("Object ID")
                or ""
            )
            obj_id = str(obj_id).strip()
            if not obj_id:
                n_skipped += 1
                continue

            geometry_data = feature.get("geometry")
            if not geometry_data:
                n_skipped += 1
                continue

            try:
                geom = _shape(geometry_data)
                if geom.is_empty or geom.geom_type not in {"Polygon", "MultiPolygon"}:
                    n_skipped += 1
                    continue
                if pixel_size_um != 1.0:
                    geom = _affinity.scale(
                        geom, xfact=pixel_size_um, yfact=pixel_size_um, origin=(0, 0)
                    )
                if not geom.is_valid:
                    geom = geom.buffer(0)
                if geom.is_empty:
                    n_skipped += 1
                    continue
            except Exception:
                n_skipped += 1
                continue

            geom_by_id[obj_id] = geom
            n_features_read += 1

    # Match to adata.obs by Object ID
    obs_ids = adata.obs[object_id_col].astype(str)
    wkt_array = _np.full(len(adata.obs), "", dtype=object)
    n_matched = 0

    for position, obj_id in enumerate(obs_ids):
        geom = geom_by_id.get(obj_id)
        if geom is not None:
            wkt_array[position] = str(geom.wkt)  # use .wkt property — avoids shapely 2.x array return
            n_matched += 1

    adata.obs[obs_key] = wkt_array
    if obsm_key is not None:
        adata.obsm[obsm_key] = wkt_array

    if verbose:
        print(
            f"load_cell_polygons: {n_features_read:,} polygons read from "
            f"{len(geojson_paths)} GeoJSON file(s), "
            f"{n_matched:,} of {len(adata.obs):,} cells matched "
            f"→ adata.obs['{obs_key}']"
        )
        if n_matched < len(adata.obs):
            n_unmatched = len(adata.obs) - n_matched
            print(
                f"  {n_unmatched:,} cells unmatched (empty string). "
                "Check that the GeoJSON was exported with Object IDs included."
            )

    return n_matched


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
