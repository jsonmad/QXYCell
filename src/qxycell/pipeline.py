"""Main QXYCell pipeline entry point."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from qxycell.classifiers import ClassifierDefinition
from qxycell.classifiers import marker_name_from_measurement_column
from qxycell.classifiers import measurement_columns_for_threshold_template
from qxycell.celltyping import apply_celltypes
from qxycell.checks import check
from qxycell.filtering import assign_core_ids_from_measurements, assign_samples
from qxycell.geojson import _classification_name
from qxycell.markers import marker_name_from_classifier_name
from qxycell.measurements import required_columns
from qxycell.paths import resolve_output_dir

CENTROID_OBS_COLUMN_RENAMES = {
    "Centroid X µm": "Xµm",
    "Centroid Y µm": "Yµm",
}


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


def _classifier_group_key(classifier: ClassifierDefinition) -> tuple[str, str]:
    return (str(classifier.name), str(classifier.measurement_column))


def _group_simple_classifiers(
    classifiers: list[ClassifierDefinition],
) -> list[list[ClassifierDefinition]]:
    grouped: dict[tuple[str, str], list[ClassifierDefinition]] = {}
    for classifier in classifiers:
        if not classifier.is_simple or classifier.measurement_column is None:
            continue
        grouped.setdefault(_classifier_group_key(classifier), []).append(classifier)
    return list(grouped.values())


def _unique_marker_names(classifier_groups: list[list[ClassifierDefinition]]) -> dict[int, str]:
    names: dict[int, str] = {}
    seen: dict[str, int] = {}
    for group_index, group in enumerate(classifier_groups):
        classifier = group[0]
        base = marker_name_from_classifier_name(classifier.name)
        count = seen.get(base, 0)
        seen[base] = count + 1
        names[group_index] = base if count == 0 else f"{base}_{count + 1}"
    return names


def _threshold_summary(classifiers: list[ClassifierDefinition]) -> str | float:
    thresholds = {float(classifier.threshold) for classifier in classifiers if classifier.threshold is not None}
    images = {classifier.image for classifier in classifiers if classifier.image}
    if len(thresholds) == 1 and not images:
        return next(iter(thresholds))
    if images:
        return "per_image"
    return "|".join(str(value) for value in sorted(thresholds))


def _build_var_dataframe(classifier_groups, marker_names, pd):
    rows = []
    index = []
    for group_index, group in enumerate(classifier_groups):
        classifier = group[0]
        marker_name = marker_names[group_index]
        index.append(marker_name)
        rows.append(
            {
                "marker_name": marker_name,
                "classifier_name": classifier.name,
                "source_measurement_column": classifier.measurement_column,
                "threshold": _threshold_summary(group),
                "threshold_source": "|".join(str(item.path) for item in group),
            }
        )
    return pd.DataFrame(rows, index=index)


def _unique_marker_names_from_measurement_columns(columns: list[str]) -> dict[int, str]:
    names: dict[int, str] = {}
    seen: dict[str, int] = {}
    for index, column in enumerate(columns):
        base = marker_name_from_measurement_column(column)
        count = seen.get(base, 0)
        seen[base] = count + 1
        names[index] = base if count == 0 else f"{base}_{count + 1}"
    return names


def _build_var_dataframe_from_measurement_columns(marker_columns, marker_names, pd):
    rows = []
    index = []
    for column_index, measurement_column in enumerate(marker_columns):
        marker_name = marker_names[column_index]
        index.append(marker_name)
        rows.append(
            {
                "marker_name": marker_name,
                "classifier_name": "",
                "source_measurement_column": measurement_column,
                "threshold": "",
                "threshold_source": "",
            }
        )
    return pd.DataFrame(rows, index=index)


def _apply_marker_thresholds(
    adata,
    classifier_groups,
    marker_names,
    image_col: str = "Image",
    marker_indices: list[int] | None = None,
) -> int:
    import numpy as _np

    n_pos_columns = 0
    for group_index, group in enumerate(classifier_groups):
        marker_name = marker_names[group_index]
        marker_index = marker_indices[group_index] if marker_indices is not None else group_index
        values = _np.asarray(adata.X[:, marker_index], dtype=float)
        pos = _np.zeros(adata.n_obs, dtype="int8")
        fallback_thresholds = [
            float(classifier.threshold)
            for classifier in group
            if classifier.threshold is not None and not classifier.image
        ]
        fallback_threshold = fallback_thresholds[-1] if fallback_thresholds else None
        if fallback_threshold is not None:
            pos[:] = (values >= fallback_threshold).astype("int8")
        for classifier in group:
            if classifier.threshold is None or not classifier.image:
                continue
            mask = adata.obs[image_col].astype(str).to_numpy() == str(classifier.image)
            pos[mask] = (values[mask] >= float(classifier.threshold)).astype("int8")
        adata.obs[f"{marker_name}_pos"] = pos
        n_pos_columns += 1
    return n_pos_columns


def _threshold_output_dir(adata, output_dir: str | Path | None = None) -> Path:
    if output_dir is not None:
        return resolve_output_dir(output_dir, adata=adata)
    metadata = getattr(adata, "uns", {}).get("qxycell", {})
    if isinstance(metadata, dict) and metadata.get("output_dir"):
        return Path(metadata["output_dir"]).expanduser().resolve()
    return resolve_output_dir(None, adata=adata)


def apply_thresholds(
    adata,
    project_dir: str | Path | None = None,
    *,
    threshold_file: str | Path | None = None,
    output_dir: str | Path | None = None,
    image_col: str = "Image",
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply marker threshold definitions to an existing AnnData object.

    This step adds ``<marker>_pos`` columns to ``adata.obs``. It is separate
    from ``qxy.run()``, which imports measurement intensities, spatial
    coordinates, and annotation metadata into AnnData.
    """

    import pandas as pd

    metadata = getattr(adata, "uns", {}).get("qxycell", {})
    if project_dir is None and isinstance(metadata, dict):
        project_dir = metadata.get("project_dir")
    if project_dir is None:
        raise ValueError(
            "project_dir is required unless adata.uns['qxycell']['project_dir'] is set."
        )
    if image_col not in adata.obs.columns:
        raise KeyError(f"image_col not found in adata.obs: {image_col}")
    if "source_measurement_column" not in adata.var.columns:
        raise KeyError("adata.var must contain 'source_measurement_column'.")

    output_path = _threshold_output_dir(adata, output_dir)
    report = check(project_dir, output_dir=output_path, threshold_file=threshold_file)
    if not report.ok:
        raise RuntimeError(
            "QXYCell check failed before thresholding. See "
            f"{output_path / 'check_report.txt'}"
        )

    simple_classifiers = [classifier for classifier in report.classifiers if classifier.is_simple]
    if not simple_classifiers:
        raise ValueError("No usable threshold definitions are available.")

    source_columns = adata.var["source_measurement_column"].astype(str).tolist()
    source_lookup = {column: index for index, column in enumerate(source_columns)}
    classifier_groups = _group_simple_classifiers(simple_classifiers)
    matched_groups = []
    marker_indices = []
    marker_names: dict[int, str] = {}
    missing_columns = []
    for group in classifier_groups:
        measurement_column = str(group[0].measurement_column)
        marker_index = source_lookup.get(measurement_column)
        if marker_index is None:
            missing_columns.append(measurement_column)
            continue
        group_index = len(matched_groups)
        matched_groups.append(group)
        marker_indices.append(marker_index)
        marker_names[group_index] = str(adata.var_names[marker_index])

    if missing_columns:
        raise ValueError(
            "Threshold definitions reference measurement columns not present in adata.var: "
            + ", ".join(sorted(dict.fromkeys(missing_columns)))
        )

    for column in ("classifier_name", "threshold", "threshold_source"):
        if column in adata.var.columns:
            adata.var[column] = adata.var[column].astype("object")

    n_pos_columns = _apply_marker_thresholds(
        adata,
        matched_groups,
        marker_names,
        image_col=image_col,
        marker_indices=marker_indices,
    )

    for group_index, group in enumerate(matched_groups):
        marker_index = marker_indices[group_index]
        adata.var.iloc[marker_index, adata.var.columns.get_loc("classifier_name")] = group[0].name
        adata.var.iloc[marker_index, adata.var.columns.get_loc("threshold")] = _threshold_summary(group)
        adata.var.iloc[marker_index, adata.var.columns.get_loc("threshold_source")] = "|".join(
            str(item.path) for item in group
        )

    summary = {
        "project_dir": str(Path(project_dir).expanduser().resolve()),
        "output_dir": str(output_path),
        "threshold_source": str(report.active_threshold_source)
        if report.active_threshold_source
        else report.active_threshold_source_kind,
        "threshold_source_kind": report.active_threshold_source_kind,
        "generated_threshold_template": str(report.generated_threshold_template)
        if report.generated_threshold_template
        else None,
        "n_threshold_definitions": len(simple_classifiers),
        "n_marker_groups": len(matched_groups),
        "n_pos_columns": n_pos_columns,
        "pos_columns": [f"{marker_names[index]}_pos" for index in range(len(matched_groups))],
    }
    adata.uns["qxycell_thresholding"] = summary
    adata.uns.setdefault("qxycell", {})["thresholding_applied"] = True
    adata.uns["qxycell"]["threshold_source"] = summary["threshold_source"]
    adata.uns["qxycell"]["threshold_source_kind"] = summary["threshold_source_kind"]
    adata.uns["qxycell"]["generated_threshold_template"] = summary["generated_threshold_template"]

    tables_dir = output_path / "run" / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(tables_dir / "thresholding_summary.csv", index=False)

    if verbose:
        print(
            "Applied marker thresholds: "
            f"{n_pos_columns} positivity columns from {summary['threshold_source']}"
        )
    return summary


threshold = apply_thresholds


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


def _load_geojson_features(
    geojson_files,
    pixel_size_um: float,
    skip_annotation_labels: set[str] | None = None,
):
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
    skip_annotation_labels = set(skip_annotation_labels or set())

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
                    if label in skip_annotation_labels:
                        continue
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


def _apply_annotations(
    adata,
    geojson_files,
    pixel_size_um: float,
    skip_annotation_labels: set[str] | None = None,
):
    try:
        from shapely.geometry import Point
    except ImportError:
        return []

    annotations_by_image = _load_geojson_features(
        geojson_files,
        pixel_size_um=pixel_size_um,
        skip_annotation_labels=skip_annotation_labels,
    )
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
    threshold_file: str | Path | None = None,
    apply_thresholds: bool = False,
    celltype_logic: str | Path | dict[str, Any] | None = None,
    verbose: bool = True,
) -> Any:
    """Run QXYCell on a manually exported QuPath project.

    The pipeline imports QuPath measurement intensity columns into ``adata.X``
    and stores required QuPath identity/spatial columns plus available core
    metadata columns in ``adata.obs``. Marker positivity calls are a separate
    step via ``qxy.threshold()`` / ``qxy.apply_thresholds()``.
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

    report = check(project_dir, output_dir=output_path, threshold_file=threshold_file)
    log(
        "Check: "
        f"{'PASS' if report.ok else 'FAIL'} "
        f"({report.n_errors} errors, {report.n_warnings} warnings)"
    )
    log(f"Measurement files: {len(report.measurement_files)}")
    log(f"Classifier definitions: {len(report.classifiers)}")
    log(f"Simple classifiers imported: {sum(1 for item in report.classifiers if item.is_simple)}")
    log(
        "Threshold source: "
        f"{report.active_threshold_source or report.active_threshold_source_kind}"
    )
    log(f"GeoJSON files: {len(report.geojson_files)}")

    if fail_on_check_error and not report.ok:
        raise RuntimeError(
            "QXYCell check failed. See "
            f"{output_path / 'check_report.txt'}"
        )

    log("Loading measurement table(s)...")
    measurements = _read_measurements(report.measurement_files, pd)
    log(f"Loaded cells: {len(measurements):,}")
    for column in required_columns():
        if column not in measurements.columns:
            raise ValueError(f"Missing required QuPath measurement column: {column}")

    marker_columns = measurement_columns_for_threshold_template(report.measurement_files)
    if not marker_columns:
        raise ValueError("No mean/median measurement columns are available for adata.X import.")
    marker_names = _unique_marker_names_from_measurement_columns(marker_columns)
    missing_marker_columns = [
        column for column in marker_columns if column not in measurements.columns
    ]
    if missing_marker_columns:
        raise ValueError(
            "Measurement columns selected for adata.X are missing: "
            + ", ".join(str(column) for column in missing_marker_columns)
        )

    log("Building marker matrix...")
    x = measurements[marker_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy()

    optional_obs_columns = [
        column
        for column in ("TMA Core", "Parent")
        if column in measurements.columns
    ]
    obs_columns = (
        list(required_columns())
        + optional_obs_columns
        + ["quxy_source_file", "quxy_source_row"]
    )
    obs = measurements.loc[:, obs_columns].copy()
    obs = obs.rename(columns=CENTROID_OBS_COLUMN_RENAMES)
    obs["quxy_cell_id"] = (
        obs["Image"].astype(str) + "::" + obs["Object ID"].astype(str)
    )
    obs.index = obs["quxy_cell_id"].astype(str)

    spatial = (
        obs[["Xµm", "Yµm"]]
        .apply(pd.to_numeric, errors="coerce")
        .to_numpy()
    )
    if np.isnan(spatial).any():
        raise ValueError("Centroid X/Y columns contain missing or non-numeric values.")

    var = _build_var_dataframe_from_measurement_columns(marker_columns, marker_names, pd)
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.var_names = list(var.index)
    adata.obsm["spatial"] = spatial

    log(f"Created AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} markers")
    log("Marker positivity columns not added by run(); call qxy.threshold(adata) next.")

    coreid_summary = None
    matching_geojson_core_labels = set(report.geojson_core_annotation_counts or {})
    if any(column in adata.obs.columns for column in ("TMA Core", "Parent")):
        coreid_summary = assign_core_ids_from_measurements(adata, verbose=False)
        log(
            "Measurement CoreID assignment: "
            f"{coreid_summary['n_assigned_cells']:,} assigned, "
            f"{coreid_summary['n_unassigned_cells']:,} unassigned"
        )
        if matching_geojson_core_labels:
            log(
                "GeoJSON core-like annotations matched measurement CoreID labels: "
                f"{len(matching_geojson_core_labels)} label(s). "
                "Measurement CoreID values were used; matching GeoJSON labels were "
                "not kept as annotation columns."
            )

    log("Mapping GeoJSON annotations...")
    annotation_conflicts = _apply_annotations(
        adata,
        report.geojson_files,
        pixel_size_um=pixel_size_um,
        skip_annotation_labels=matching_geojson_core_labels if coreid_summary else set(),
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
        log(
            "TMA core features detected: "
            f"{n_tma_core_features}. "
            "Run qxy.assign_tma_cores(...) explicitly if geometry-based TMA "
            "assignment is needed."
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

    thresholding_summary = None
    if apply_thresholds:
        log("Applying marker thresholds...")
        thresholding_summary = globals()["apply_thresholds"](
            adata,
            project_dir=project_dir,
            threshold_file=threshold_file,
            output_dir=output_path,
            verbose=verbose,
        )
        log(f"Thresholding: {thresholding_summary['n_pos_columns']} positivity columns")

    celltyping_summary = None
    if celltype_logic is not None:
        if not any(str(column).endswith("_pos") for column in adata.obs.columns):
            raise ValueError(
                "celltype_logic requires marker positivity columns. "
                "Call qxy.threshold(adata, ...) first or pass apply_thresholds=True."
            )
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
        "threshold_source": str(report.active_threshold_source)
        if report.active_threshold_source
        else report.active_threshold_source_kind,
        "threshold_source_kind": report.active_threshold_source_kind,
        "generated_threshold_template": str(report.generated_threshold_template)
        if report.generated_threshold_template
        else None,
        "thresholding_applied": bool(thresholding_summary is not None),
        "thresholding": thresholding_summary,
        "n_geojson_files": int(len(report.geojson_files)),
        "n_tma_core_features": int(n_tma_core_features),
        "n_measurement_core_labels": int(report.n_measurement_core_labels),
        "n_measurement_core_cells": int(report.n_measurement_core_cells),
        "n_geojson_core_annotation_features": int(
            report.n_geojson_core_annotation_features
        ),
        "measurement_core_counts": dict(report.measurement_core_counts or {}),
        "geojson_core_annotation_counts": dict(
            report.geojson_core_annotation_counts or {}
        ),
        "measurement_core_assignment": coreid_summary,
        "ignored_geojson_core_annotation_labels": sorted(matching_geojson_core_labels)
        if coreid_summary
        else [],
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
