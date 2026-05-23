"""Main QUXYCell pipeline entry point."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from quxycell.classifiers import ClassifierDefinition
from quxycell.checks import check
from quxycell.geojson import _classification_name
from quxycell.markers import marker_name_from_classifier_name
from quxycell.measurements import required_columns


@dataclass(frozen=True)
class RunResult:
    """Result returned by :func:`run`."""

    adata: Any
    output_dir: Path
    h5ad_path: Path
    check_report_path: Path


def _import_runtime_dependencies():
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            "QUXYCell run() requires the package runtime dependencies. "
            "Install with `pip install -e .` from this repo, or `pip install quxycell` "
            "once the package is published."
        ) from exc
    return ad, np, pd


def _read_measurements(measurement_files, pd):
    frames = []
    source_columns = []
    for measurement_file in measurement_files:
        frame = pd.read_csv(measurement_file.path, sep=measurement_file.delimiter, low_memory=False)
        frame["quxy_source_file"] = str(measurement_file.path)
        frame["quxy_source_row"] = range(1, len(frame) + 1)
        frames.append(frame)
        source_columns.extend(measurement_file.columns)
    if not frames:
        raise ValueError("No measurement files available to load.")
    return pd.concat(frames, axis=0, ignore_index=True, sort=False)


def _unique_marker_names(classifiers: list[ClassifierDefinition]) -> dict[Path, str]:
    names: dict[Path, str] = {}
    seen: dict[str, int] = {}
    for classifier in classifiers:
        if not classifier.is_simple:
            continue
        base = marker_name_from_classifier_name(classifier.name)
        count = seen.get(base, 0)
        seen[base] = count + 1
        names[classifier.path] = base if count == 0 else f"{base}_{count + 1}"
    return names


def _build_var_dataframe(simple_classifiers, marker_names, pd):
    rows = []
    index = []
    for classifier in simple_classifiers:
        marker_name = marker_names[classifier.path]
        index.append(marker_name)
        rows.append(
            {
                "marker_name": marker_name,
                "classifier_name": classifier.name,
                "source_measurement_column": classifier.measurement_column,
                "threshold": classifier.threshold,
                "classifier_json": str(classifier.path),
            }
        )
    return pd.DataFrame(rows, index=index)


def _safe_obs_column_name(prefix: str, label: str) -> str:
    safe = "".join(ch if ch.isalnum() else "_" for ch in str(label).strip())
    safe = "_".join(part for part in safe.split("_") if part)
    return f"{prefix}__{safe or 'Unclassified'}"


def _geojson_image_key(path: Path) -> str:
    stem = path.stem
    return stem[:-4] if stem.endswith(".ome") else stem


def _image_key(value: str) -> str:
    path = Path(str(value))
    stem = path.stem
    return stem[:-4] if stem.endswith(".ome") else stem


def _load_annotation_features(geojson_files, pixel_size_um: float):
    try:
        import json
        from shapely import affinity
        from shapely.geometry import shape
        from shapely.prepared import prep
    except ImportError as exc:
        raise ImportError(
            "GeoJSON annotation mapping requires shapely. Install QUXYCell runtime dependencies."
        ) from exc

    by_image: dict[str, list[dict[str, Any]]] = {}
    for geojson_file in geojson_files:
        if not geojson_file.readable:
            continue
        data = json.loads(geojson_file.path.read_text(errors="replace"))
        features = data.get("features", []) if isinstance(data, dict) else []
        image_key = _geojson_image_key(geojson_file.path)
        for feature in features:
            properties = feature.get("properties", {}) if isinstance(feature, dict) else {}
            geometry_data = feature.get("geometry") if isinstance(feature, dict) else None
            if not geometry_data:
                continue
            geometry = shape(geometry_data)
            if geometry.is_empty or geometry.geom_type not in {"Polygon", "MultiPolygon"}:
                continue
            if pixel_size_um != 1.0:
                geometry = affinity.scale(
                    geometry,
                    xfact=pixel_size_um,
                    yfact=pixel_size_um,
                    origin=(0, 0),
                )
            if not geometry.is_valid:
                geometry = geometry.buffer(0)
            if geometry.is_empty:
                continue
            label = _classification_name(properties)
            by_image.setdefault(image_key, []).append(
                {
                    "label": label,
                    "column": _safe_obs_column_name("annotation", label),
                    "prepared": prep(geometry),
                    "source": str(geojson_file.path),
                }
            )
    return by_image


def _apply_annotations(adata, geojson_files, pixel_size_um: float):
    try:
        from shapely.geometry import Point
    except ImportError:
        return []

    features_by_image = _load_annotation_features(geojson_files, pixel_size_um=pixel_size_um)
    if not features_by_image:
        return []

    obs = adata.obs
    annotation_columns = sorted(
        {feature["column"] for features in features_by_image.values() for feature in features}
    )
    for column in annotation_columns:
        obs[column] = False

    spatial = adata.obsm["spatial"]
    image_keys = obs["Image"].astype(str).map(_image_key)
    conflict_rows = []
    for image_key, features in features_by_image.items():
        indices = list(obs.index[image_keys == image_key])
        if not indices:
            continue
        positions = [obs.index.get_loc(index) for index in indices]
        for index, position in zip(indices, positions, strict=True):
            point = Point(float(spatial[position, 0]), float(spatial[position, 1]))
            label_hits: dict[str, int] = {}
            for feature in features:
                if feature["prepared"].contains(point):
                    obs.at[index, feature["column"]] = True
                    label_hits[feature["column"]] = label_hits.get(feature["column"], 0) + 1
            for column, count in label_hits.items():
                if count > 1:
                    conflict_rows.append(
                        {
                            "Image": obs.at[index, "Image"],
                            "Object ID": obs.at[index, "Object ID"],
                            "annotation_column": column,
                            "hit_count": count,
                        }
                    )
    return conflict_rows


def run(
    project_dir: str | Path,
    output_dir: str | Path = "quxycell_output",
    *,
    fail_on_check_error: bool = True,
    pixel_size_um: float = 0.28,
) -> RunResult:
    """Run QUXYCell on a manually exported QuPath project.

    The v1 pipeline imports one ``adata.X`` column per usable simple measurement
    classifier JSON, stores only required QuPath identity/spatial columns in
    ``adata.obs``, and adds classifier positivity calls plus GeoJSON annotation
    columns when possible.
    """

    ad, np, pd = _import_runtime_dependencies()

    output_path = Path(output_dir).expanduser().resolve()
    report = check(project_dir, output_dir=output_path)
    if fail_on_check_error and not report.ok:
        raise RuntimeError(
            "QUXYCell check failed. See "
            f"{output_path / 'check_report.txt'}"
        )

    simple_classifiers = [classifier for classifier in report.classifiers if classifier.is_simple]
    if not simple_classifiers:
        raise ValueError("No simple classifier JSON files are available for adata.X import.")

    measurements = _read_measurements(report.measurement_files, pd)
    for column in required_columns():
        if column not in measurements.columns:
            raise ValueError(f"Missing required QuPath measurement column: {column}")

    marker_names = _unique_marker_names(simple_classifiers)
    marker_columns = [classifier.measurement_column for classifier in simple_classifiers]
    missing_marker_columns = [column for column in marker_columns if column not in measurements.columns]
    if missing_marker_columns:
        raise ValueError(
            "Classifier-referenced measurement columns are missing: "
            + ", ".join(str(column) for column in missing_marker_columns)
        )

    x = measurements[marker_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy()

    obs_columns = list(required_columns()) + ["quxy_source_file", "quxy_source_row"]
    obs = measurements.loc[:, obs_columns].copy()
    obs["quxy_cell_id"] = (
        obs["Image"].astype(str) + "::" + obs["Object ID"].astype(str)
    )
    obs.index = obs["quxy_cell_id"].astype(str)

    spatial = (
        obs[["Centroid X µm", "Centroid Y µm"]]
        .apply(pd.to_numeric, errors="coerce")
        .to_numpy()
    )
    if np.isnan(spatial).any():
        raise ValueError("Centroid X/Y columns contain missing or non-numeric values.")

    var = _build_var_dataframe(simple_classifiers, marker_names, pd)
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.var_names = list(var.index)
    adata.obsm["spatial"] = spatial

    for classifier_index, classifier in enumerate(simple_classifiers):
        marker_name = marker_names[classifier.path]
        adata.obs[f"{marker_name}_pos"] = (
            adata.X[:, classifier_index] >= float(classifier.threshold)
        ).astype("int8")

    annotation_conflicts = _apply_annotations(
        adata,
        report.geojson_files,
        pixel_size_um=pixel_size_um,
    )
    adata.uns["quxycell"] = {
        "project_dir": str(Path(project_dir).expanduser().resolve()),
        "output_dir": str(output_path),
        "created": datetime.now().isoformat(timespec="seconds"),
        "pixel_size_um": pixel_size_um,
        "check_report_txt": str(output_path / "check_report.txt"),
        "check_report_json": str(output_path / "check_report.json"),
        "check_ok": bool(report.ok),
        "check_n_errors": int(report.n_errors),
        "check_n_warnings": int(report.n_warnings),
        "n_measurement_files": int(len(report.measurement_files)),
        "n_classifiers": int(len(report.classifiers)),
        "n_simple_classifiers": int(sum(1 for item in report.classifiers if item.is_simple)),
        "n_geojson_files": int(len(report.geojson_files)),
        "n_annotation_conflicts": int(len(annotation_conflicts)),
    }

    run_dir = output_path / "run"
    h5ad_dir = run_dir / "h5ad"
    h5ad_dir.mkdir(parents=True, exist_ok=True)
    h5ad_path = h5ad_dir / "quxycell.h5ad"
    adata.write_h5ad(h5ad_path)

    tables_dir = run_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    adata.obs.to_csv(tables_dir / "cells_obs.csv")
    adata.var.to_csv(tables_dir / "markers_var.csv")
    pd.DataFrame(annotation_conflicts).to_csv(
        tables_dir / "annotation_conflicts.csv",
        index=False,
    )

    return RunResult(
        adata=adata,
        output_dir=run_dir,
        h5ad_path=h5ad_path,
        check_report_path=output_path / "check_report.txt",
    )
