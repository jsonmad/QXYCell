"""Main QXYCell pipeline entry point."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from qxycell.classifiers import ClassifierDefinition
from qxycell.celltyping import apply_celltypes
from qxycell.checks import check
from qxycell.filtering import assign_samples
from qxycell.geojson import _classification_name
from qxycell.markers import marker_name_from_classifier_name
from qxycell.measurements import required_columns
from qxycell.paths import resolve_output_dir


def _import_runtime_dependencies():
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            "QXYCell run() requires the package runtime dependencies. "
            "Install with `pip install -e .` from this repo, or `pip install qxycell` "
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


def _is_sample_annotation_label(label: str) -> bool:
    return "sample" in str(label).lower()


def _annotation_labels(properties: dict[str, Any]) -> list[str]:
    labels = [_classification_name(properties)]
    name = properties.get("name")
    if name:
        label = str(name)
        if label and label.lower() not in {"none", "null"} and label not in labels:
            labels.append(label)
    return labels


def _geojson_image_key(path: Path) -> str:
    stem = path.stem
    return stem[:-4] if stem.endswith(".ome") else stem


def _image_key(value: str) -> str:
    path = Path(str(value))
    stem = path.stem
    return stem[:-4] if stem.endswith(".ome") else stem


def _load_geojson_features(geojson_files, pixel_size_um: float):
    """Load annotation features from GeoJSON files, keyed by image.

    Returns:
      {image_key: [{"label", "column", "prepared", "source"}, ...]}
    """
    try:
        import json
        from shapely import affinity
        from shapely.geometry import shape
        from shapely.prepared import prep
    except ImportError as exc:
        raise ImportError(
            "GeoJSON annotation mapping requires shapely. Install QXYCell runtime dependencies."
        ) from exc

    def _prepare(feature_dict) -> Any:
        geometry_data = feature_dict.get("geometry")
        if not geometry_data:
            return None
        geometry = shape(geometry_data)
        if geometry.is_empty or geometry.geom_type not in {"Polygon", "MultiPolygon"}:
            return None
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
            return None
        return prep(geometry)

    annotations_by_image: dict[str, list[dict[str, Any]]] = {}

    for geojson_file in geojson_files:
        if not geojson_file.readable:
            continue
        data = json.loads(geojson_file.path.read_text(errors="replace"))
        features = data.get("features", []) if isinstance(data, dict) else []
        image_key = _geojson_image_key(geojson_file.path)

        for feature in features:
            if not isinstance(feature, dict):
                continue
            properties = feature.get("properties") or {}
            if not isinstance(properties, dict):
                properties = {}
            object_type = str(properties.get("objectType") or "").lower()

            if object_type == "annotation":
                prepared = _prepare(feature)
                if prepared is None:
                    continue
                for label in _annotation_labels(properties):
                    prefix = "sample_annotation" if _is_sample_annotation_label(label) else "annotation"
                    annotations_by_image.setdefault(image_key, []).append(
                        {
                            "label": label,
                            "column": _safe_obs_column_name(prefix, label),
                            "is_sample": _is_sample_annotation_label(label),
                            "prepared": prepared,
                            "source": str(geojson_file.path),
                        }
                    )

    return annotations_by_image


def _apply_cell_polygons(adata, geojson_files, pixel_size_um: float) -> int:
    """Match cell detection polygons from GeoJSON to adata.obs rows via Object ID.

    Stores WKT geometry strings in ``adata.obs["cell_polygon_wkt"]`` (empty string for
    unmatched cells). Returns the number of cells successfully matched.
    """
    try:
        import json as _json
        import numpy as _np
        from shapely import affinity
        from shapely.geometry import shape
    except ImportError:
        return 0

    # QuPath Object IDs are globally unique UUIDs. Match by Object ID alone so
    # standard exports such as "image-cells.geojson" do not depend on filename
    # heuristics to recover the measurement-table Image value.
    cell_geoms: dict[str, Any] = {}

    for geojson_file in geojson_files:
        if not geojson_file.readable:
            continue
        data = _json.loads(geojson_file.path.read_text(errors="replace"))
        features = data.get("features", []) if isinstance(data, dict) else []

        for feature in features:
            if not isinstance(feature, dict):
                continue
            properties = feature.get("properties") or {}
            if not isinstance(properties, dict):
                properties = {}
            object_type = str(properties.get("objectType") or "").lower()
            if object_type != "cell":
                continue
            geometry_data = feature.get("geometry")
            if not geometry_data:
                continue
            try:
                geom = shape(geometry_data)
                if geom.is_empty or geom.geom_type not in {"Polygon", "MultiPolygon"}:
                    continue
                if pixel_size_um != 1.0:
                    geom = affinity.scale(geom, xfact=pixel_size_um, yfact=pixel_size_um, origin=(0, 0))
                if not geom.is_valid:
                    geom = geom.buffer(0)
                if geom.is_empty:
                    continue
            except Exception:
                continue

            # QuPath stores Object ID as the top-level GeoJSON feature "id",
            # not inside "properties". Fall back to property keys for others.
            obj_id = (
                feature.get("id")
                or properties.get("objectID")
                or properties.get("id")
                or properties.get("Object ID")
                or ""
            )
            obj_id = str(obj_id).strip()
            if obj_id:
                cell_geoms[obj_id] = geom

    if not cell_geoms:
        return 0

    obs = adata.obs
    object_ids = obs["Object ID"].astype(str)

    wkt_array = _np.full(len(obs), "", dtype=object)
    n_matched = 0

    for position, obj_id in enumerate(object_ids):
        geom = cell_geoms.get(obj_id)
        if geom is not None:
            wkt_array[position] = str(geom.wkt)
            n_matched += 1

    adata.obs["cell_polygon_wkt"] = wkt_array
    return n_matched


def _apply_annotations(adata, geojson_files, pixel_size_um: float):
    try:
        from shapely.geometry import Point
    except ImportError:
        return []

    annotations_by_image = _load_geojson_features(geojson_files, pixel_size_um=pixel_size_um)
    if not annotations_by_image:
        return []

    obs = adata.obs

    # --- Annotation boolean columns ---
    annotation_columns = sorted(
        {feature["column"] for features in annotations_by_image.values() for feature in features}
    )
    annotation_label_map = {
        feature["column"]: feature["label"]
        for features in annotations_by_image.values()
        for feature in features
    }
    adata.uns["qxycell_annotation_labels"] = annotation_label_map
    for column in annotation_columns:
        obs[column] = False

    spatial = adata.obsm["spatial"]
    image_keys = obs["Image"].astype(str).map(_image_key)

    for image_key, annotation_features in annotations_by_image.items():
        indices = list(obs.index[image_keys == image_key])
        if not indices:
            continue
        positions = [obs.index.get_loc(index) for index in indices]

        for index, position in zip(indices, positions, strict=True):
            point = Point(float(spatial[position, 0]), float(spatial[position, 1]))

            # Annotation assignment
            for feature in annotation_features:
                if feature["prepared"].contains(point):
                    obs.at[index, feature["column"]] = True

    return []


def run(
    project_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    fail_on_check_error: bool = True,
    pixel_size_um: float = 0.28,
    celltype_logic: str | Path | dict[str, Any] | None = None,
    verbose: bool = True,
) -> Any:
    """Run QXYCell on a manually exported QuPath project.

    The v1 pipeline imports one ``adata.X`` column per usable simple measurement
    classifier JSON, stores only required QuPath identity/spatial columns in
    ``adata.obs``, and adds classifier positivity calls plus GeoJSON annotation
    columns when possible.
    """

    ad, np, pd = _import_runtime_dependencies()

    output_path = resolve_output_dir(output_dir)
    log_lines: list[str] = []

    def log(message: str) -> None:
        log_lines.append(message)
        if verbose:
            print(message)

    log("QXYCell run started")
    log(f"Project: {Path(project_dir).expanduser().resolve()}")
    log(f"Output: {output_path}")

    report = check(project_dir, output_dir=output_path)
    log(
        "Check: "
        f"{'PASS' if report.ok else 'FAIL'} "
        f"({report.n_errors} errors, {report.n_warnings} warnings)"
    )
    log(f"Measurement files: {len(report.measurement_files)}")
    log(f"Classifier JSON files: {len(report.classifiers)}")
    log(f"Simple classifiers imported: {sum(1 for item in report.classifiers if item.is_simple)}")
    log(f"GeoJSON files: {len(report.geojson_files)}")

    if fail_on_check_error and not report.ok:
        raise RuntimeError(
            "QXYCell check failed. See "
            f"{output_path / 'check_report.txt'}"
        )

    simple_classifiers = [classifier for classifier in report.classifiers if classifier.is_simple]
    if not simple_classifiers:
        raise ValueError("No simple classifier JSON files are available for adata.X import.")

    log("Loading measurement table(s)...")
    measurements = _read_measurements(report.measurement_files, pd)
    log(f"Loaded cells: {len(measurements):,}")
    for column in required_columns():
        if column not in measurements.columns:
            raise ValueError(f"Missing required QuPath measurement column: {column}")

    marker_names = _unique_marker_names(simple_classifiers)
    marker_columns = [classifier.measurement_column for classifier in simple_classifiers]
    missing_marker_columns = [
        column for column in marker_columns if column not in measurements.columns
    ]
    if missing_marker_columns:
        raise ValueError(
            "Classifier-referenced measurement columns are missing: "
            + ", ".join(str(column) for column in missing_marker_columns)
        )

    log("Building marker matrix...")
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
    log(f"Created AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} markers")
    log(f"Added marker positivity columns: {len(simple_classifiers)}")

    log("Mapping GeoJSON annotations...")
    annotation_conflicts = _apply_annotations(
        adata,
        report.geojson_files,
        pixel_size_um=pixel_size_um,
    )
    annotation_cols = [
        column for column in adata.obs.columns if str(column).startswith("annotation__")
    ]
    sample_annotation_cols = [
        column for column in adata.obs.columns if str(column).startswith("sample_annotation__")
    ]
    annotation_label_map = adata.uns.get("qxycell_annotation_labels", {})
    log(f"Annotation columns: {len(annotation_cols)}")
    log(f"Annotation conflicts: {len(annotation_conflicts)}")
    sample_summary = None
    if sample_annotation_cols:
        sample_summary = assign_samples(
            adata,
            annotation_prefix="sample_annotation__",
            sample_col="Sample",
            verbose=verbose,
        )
        sample_summary_for_uns = adata.uns.get("qxycell_sample_annotations", sample_summary)
        adata.obs.drop(columns=sample_annotation_cols, inplace=True)
        annotation_label_map = {
            column: label
            for column, label in annotation_label_map.items()
            if not str(column).startswith("sample_annotation__")
        }
        adata.uns["qxycell_annotation_labels"] = annotation_label_map
        log(
            "Sample annotations: "
            f"{len(sample_annotation_cols)} labels, "
            f"{sample_summary['n_assigned_cells']:,} assigned, "
            f"{sample_summary['n_conflicting_cells']:,} ambiguous"
        )

    n_tma_core_features = sum(
        count
        for geojson_file in report.geojson_files
        for object_type, count in geojson_file.object_type_counts.items()
        if str(object_type).lower() == "tmacore"
    )
    tma_summary = None
    if n_tma_core_features:
        from qxycell.tma import assign_tma_cores

        log(
            "TMA core features detected: "
            f"{n_tma_core_features}. "
            "Assigning TMA cores with qxy.assign_tma_cores(...)."
        )
        tma_summary = assign_tma_cores(
            adata,
            Path(project_dir).expanduser().resolve(),
            pixel_size_um=pixel_size_um,
            output_dir=output_path,
            verbose=verbose,
        )
        log(
            "TMA assignment: "
            f"{tma_summary['n_cores']} cores, "
            f"{tma_summary['n_assigned_cells']:,}/{tma_summary['n_cells']:,} cells assigned"
        )

    log("Mapping cell segmentation polygons...")
    n_matched = _apply_cell_polygons(adata, report.geojson_files, pixel_size_um=pixel_size_um)
    if n_matched > 0:
        log(
            "Cell polygons matched: "
            f"{n_matched:,} of {adata.n_obs:,} cells → "
            "adata.obs['cell_polygon_wkt']"
        )
    else:
        log("Cell segmentation GeoJSON not found or no Object ID matches — skipping.")

    # Per-image cell counts.
    image_counts = adata.obs["Image"].value_counts().sort_index()
    log("")
    log("Cells per image:")
    for image_name, count in image_counts.items():
        log(f"  {image_name}: {count:,}")

    # Ignore annotation summary.
    ignore_col = "annotation__Ignore"
    if ignore_col in adata.obs.columns:
        n_ignore = int(adata.obs[ignore_col].sum())
        log(f"Cells inside Ignore region(s): {n_ignore:,} of {adata.n_obs:,} "
            f"({100 * n_ignore / adata.n_obs:.1f}%) — remove with qxy.remove_ignore(adata)")
    else:
        log("No Ignore annotations found.")

    log("")

    celltyping_summary = None
    if celltype_logic is not None:
        log("Applying cell type logic...")
        celltyping_summary = apply_celltypes(
            adata,
            celltype_logic,
            celltype_dir=output_path / "celltype",
        )
        log(
            "Celltyping: "
            f"{celltyping_summary['n_rules']} rules, "
            f"{celltyping_summary['unknown_count']:,} Unknown cells"
        )

    run_dir = output_path / "run"
    h5ad_dir = run_dir / "h5ad"
    h5ad_dir.mkdir(parents=True, exist_ok=True)
    _folder_name = output_path.name  # e.g. qxy_outputs_260527-2029
    _ts = _folder_name.removeprefix("qxy_outputs_") if _folder_name.startswith("qxy_outputs_") else ""
    _h5ad_stem = f"qxycell_{_ts}" if _ts else "qxycell"
    h5ad_path = h5ad_dir / f"{_h5ad_stem}.h5ad"
    tables_dir = run_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    adata.uns["qxycell"] = {
        "project_dir": str(Path(project_dir).expanduser().resolve()),
        "output_dir": str(output_path),
        "run_dir": str(run_dir),
        "h5ad_path": str(h5ad_path),
        "tables_dir": str(tables_dir),
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
        "n_tma_core_features": int(n_tma_core_features),
        "tma_assignment": tma_summary,
        "n_annotation_conflicts": int(len(annotation_conflicts)),
        "annotation_labels": dict(annotation_label_map),
        "sample_assignment": sample_summary_for_uns if sample_summary is not None else None,
        "celltyping_applied": bool(celltyping_summary is not None),
    }
    log(f"Writing H5AD: {h5ad_path}")
    adata.write_h5ad(h5ad_path)

    log(f"Writing tables: {tables_dir}")
    adata.obs.to_csv(tables_dir / "cells_obs.csv")
    adata.var.to_csv(tables_dir / "markers_var.csv")
    pd.DataFrame(annotation_conflicts).to_csv(
        tables_dir / "annotation_conflicts.csv",
        index=False,
    )
    if celltyping_summary is not None:
        counts = adata.obs["celltype"].value_counts(dropna=False).rename_axis("celltype")
        counts.reset_index(name="cell_count").to_csv(
            tables_dir / "celltype_counts.csv",
            index=False,
        )
    log("QXYCell run complete")
    (run_dir / "run.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    return adata
