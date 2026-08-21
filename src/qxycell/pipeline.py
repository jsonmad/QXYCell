"""Main QXYCell pipeline entry point."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from qxycell.classifiers import ClassifierDefinition
from qxycell.classifiers import classifier_threshold_conflicts
from qxycell.classifiers import discover_classifier_files
from qxycell.classifiers import marker_name_from_measurement_column
from qxycell.classifiers import measurement_columns_for_threshold_template
from qxycell.classifiers import parse_classifiers
from qxycell.classifiers import parse_threshold_files
from qxycell.classifiers import unresolved_threshold_conflicts
from qxycell.checks import inspect_project, write_classifier_threshold_table
from qxycell.filtering import assign_core_ids_from_measurements, assign_samples
from qxycell.geojson import (
    _classification_name,
    discover_geojson_files,
    summarize_geojson_files,
    validate_geojson_files,
)
from qxycell.markers import marker_name_from_classifier, marker_name_from_classifier_name
from qxycell.measurements import (
    MEASUREMENT_COLUMN_ALIASES,
    MEASUREMENT_TEXT_ENCODING,
    discover_measurement_files,
    required_columns,
    summarize_measurement_file,
    validate_measurement_files,
)
from qxycell.paths import resolve_output_dir
from qxycell.stage_state import complete_stage, prepare_stage

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
            "QXYCell requires the package runtime dependencies. "
            "Install with `pip install .` from this repo, or `pip install qxycell` "
            "once the package is published."
        ) from exc
    return ad, np, pd


def _validate_pixel_size_um(value: object) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("pixel_size_um must be a positive finite number") from exc
    if not math.isfinite(result) or result <= 0:
        raise ValueError("pixel_size_um must be a positive finite number")
    return result


def _read_measurements(measurement_files, pd):
    frames = []
    source_columns = []
    for measurement_file in measurement_files:
        frame = pd.read_csv(
            measurement_file.path,
            sep=measurement_file.delimiter,
            low_memory=False,
            encoding=MEASUREMENT_TEXT_ENCODING,
        )
        for alias, canonical in MEASUREMENT_COLUMN_ALIASES.items():
            if alias not in frame.columns:
                continue
            if canonical in frame.columns:
                frame = frame.drop(columns=alias)
            else:
                frame = frame.rename(columns={alias: canonical})
        frame["quxy_source_file"] = str(measurement_file.path)
        frame["quxy_source_row"] = range(1, len(frame) + 1)
        frames.append(frame)
        source_columns.extend(measurement_file.columns)
    if not frames:
        raise ValueError("No measurement files available to load.")
    return pd.concat(frames, axis=0, ignore_index=True, sort=False)


def _classifier_group_key(classifier: ClassifierDefinition) -> str:
    return str(classifier.measurement_column)


def _group_simple_classifiers(
    classifiers: list[ClassifierDefinition],
) -> list[list[ClassifierDefinition]]:
    grouped: dict[str, list[ClassifierDefinition]] = {}
    for classifier in classifiers:
        if not classifier.is_simple or classifier.measurement_column is None:
            continue
        grouped.setdefault(_classifier_group_key(classifier), []).append(classifier)
    return list(grouped.values())


def _unique_marker_names(classifier_groups: list[list[ClassifierDefinition]]) -> dict[int, str]:
    bases = [
        marker_name_from_classifier(group[0].name, group[0].measurement_column)
        for group in classifier_groups
    ]
    return _allocate_unique_names(bases)


def _allocate_unique_names(bases: list[str]) -> dict[int, str]:
    """Allocate deterministic names without colliding with earlier outputs."""

    names: dict[int, str] = {}
    used: set[str] = set()
    for index, base in enumerate(bases):
        candidate = base
        suffix = 2
        while candidate in used:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used.add(candidate)
        names[index] = candidate
    return names


def _threshold_summary(classifiers: list[ClassifierDefinition]) -> str | float:
    thresholds = {float(classifier.threshold) for classifier in classifiers if classifier.threshold is not None}
    images = {classifier.image for classifier in classifiers if classifier.image}
    if len(thresholds) == 1 and not images:
        return next(iter(thresholds))
    if images:
        return "per_image"
    return "|".join(str(value) for value in sorted(thresholds))


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
                "threshold_marker_name": "",
                "positivity_column": "",
            }
        )
    return pd.DataFrame(rows, index=index)


def _marker_name_lookup_by_measurement_column(
    classifiers: list[ClassifierDefinition],
) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for classifier in classifiers:
        if not classifier.is_simple or classifier.measurement_column is None:
            continue
        lookup.setdefault(
            str(classifier.measurement_column),
            marker_name_from_classifier_name(classifier.name),
        )
    return lookup


def _unique_marker_names_for_measurement_columns(
    columns: list[str],
    classifiers: list[ClassifierDefinition],
) -> dict[int, str]:
    classifier_marker_lookup = _marker_name_lookup_by_measurement_column(classifiers)
    bases = [
        classifier_marker_lookup.get(
            str(column),
            marker_name_from_measurement_column(column),
        )
        for column in columns
    ]
    return _allocate_unique_names(bases)


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
    from ``qxy.import_cells()``, which imports measurement intensities and
    spatial coordinates into AnnData.
    """

    metadata = getattr(adata, "uns", {}).get("qxycell", {})
    if project_dir is None and isinstance(metadata, dict):
        project_dir = metadata.get("project_dir")
    if project_dir is None:
        raise ValueError(
            "project_dir is required unless adata.uns['qxycell']['project_dir'] is set."
        )
    output_path = _threshold_output_dir(adata, output_dir)
    report = inspect_project(project_dir, threshold_file=threshold_file, output_dir=output_path)
    if not report.ok:
        raise RuntimeError(
            "QXYCell project validation failed before thresholding. "
            "Run qxy.check(...) for a detailed report."
        )

    if report.active_threshold_source_kind == "object_classifiers":
        conflicts = classifier_threshold_conflicts(report.classifiers)
        if conflicts:
            measurements = ", ".join(
                str(item["measurement_column"]) for item in conflicts
            )
            raise ValueError(
                "Conflicting object-classifier thresholds must be resolved in a "
                "threshold table before application. Generate a table with "
                "qxy.generate_threshold_table(...), fill every image value for: "
                + measurements
            )
    if report.active_threshold_source is not None:
        unresolved = unresolved_threshold_conflicts(report.active_threshold_source)
        if unresolved:
            details = "; ".join(
                f"{item['marker'] or item['measurement_column']}: "
                + ", ".join(item["missing_images"])
                for item in unresolved
            )
            raise ValueError(
                "Threshold table contains unresolved classifier conflicts. "
                "Enter a threshold for every listed image before applying: "
                + details
            )

    return _apply_threshold_definitions(
        adata,
        report.classifiers,
        project_dir=project_dir,
        output_path=output_path,
        threshold_source=(
            str(report.active_threshold_source)
            if report.active_threshold_source
            else report.active_threshold_source_kind
        ),
        threshold_source_kind=report.active_threshold_source_kind,
        generated_threshold_template=(
            str(report.generated_threshold_template)
            if report.generated_threshold_template
            else None
        ),
        image_col=image_col,
        verbose=verbose,
    )


def _apply_threshold_definitions(
    adata,
    classifiers: list[ClassifierDefinition],
    *,
    project_dir: str | Path,
    output_path: Path,
    threshold_source: str,
    threshold_source_kind: str,
    generated_threshold_template: str | None = None,
    threshold_table_writer: Any | None = None,
    image_col: str = "Image",
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply already-selected definitions and replace prior threshold outputs."""

    import pandas as pd

    if image_col not in adata.obs.columns:
        raise KeyError(f"image_col not found in adata.obs: {image_col}")
    if "source_measurement_column" not in adata.var.columns:
        raise KeyError("adata.var must contain 'source_measurement_column'.")
    simple_classifiers = [classifier for classifier in classifiers if classifier.is_simple]
    if not simple_classifiers:
        raise ValueError("No usable threshold definitions are available.")

    source_columns = adata.var["source_measurement_column"].astype(str).tolist()
    source_lookup = {column: index for index, column in enumerate(source_columns)}
    classifier_groups = _group_simple_classifiers(simple_classifiers)
    matched_groups = []
    marker_indices = []
    missing_columns = []
    for group in classifier_groups:
        measurement_column = str(group[0].measurement_column)
        marker_index = source_lookup.get(measurement_column)
        if marker_index is None:
            missing_columns.append(measurement_column)
            continue
        matched_groups.append(group)
        marker_indices.append(marker_index)
    if missing_columns:
        raise ValueError(
            "Threshold definitions reference measurement columns not present in adata.var: "
            + ", ".join(sorted(dict.fromkeys(missing_columns)))
        )

    matched = sorted(
        zip(marker_indices, matched_groups, strict=True),
        key=lambda item: item[0],
    )
    marker_indices = [marker_index for marker_index, _group in matched]
    matched_groups = [group for _marker_index, group in matched]
    marker_names = _unique_marker_names(matched_groups)
    active_pos_columns = [
        f"{marker_names[index]}_pos" for index in range(len(matched_groups))
    ]
    metadata = adata.uns.setdefault("qxycell", {})
    stages = metadata.setdefault("stages", {})
    previous_thresholding = adata.uns.get("qxycell_thresholding", {})
    if "thresholds" not in stages and isinstance(previous_thresholding, dict):
        stages["thresholds"] = {
            "status": "complete",
            "columns": list(previous_thresholding.get("pos_columns", [])),
            "files": [],
        }
    previous_celltyping = adata.uns.get("qxycell_celltyping", {})
    if "celltypes" not in stages and isinstance(previous_celltyping, dict):
        celltype_columns = [
            previous_celltyping.get("celltype_column"),
            *previous_celltyping.get("feature_columns", []),
            *previous_celltyping.get("derived_feature_columns", []),
        ]
        stages["celltypes"] = {
            "status": "complete",
            "columns": [column for column in celltype_columns if column],
            "files": [],
        }
    replacement = prepare_stage(adata, "thresholds")
    if threshold_table_writer is not None:
        generated_threshold_template = str(
            Path(threshold_table_writer()).expanduser().resolve()
        )

    for column in (
        "classifier_name",
        "threshold",
        "threshold_source",
        "threshold_marker_name",
        "positivity_column",
    ):
        if column in adata.var.columns:
            adata.var[column] = ""
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
        adata.var.iloc[marker_index, adata.var.columns.get_loc("threshold")] = str(
            _threshold_summary(group)
        )
        adata.var.iloc[marker_index, adata.var.columns.get_loc("threshold_source")] = "|".join(
            str(item.path) for item in group
        )
        adata.var.iloc[
            marker_index, adata.var.columns.get_loc("threshold_marker_name")
        ] = marker_names[group_index]
        adata.var.iloc[
            marker_index, adata.var.columns.get_loc("positivity_column")
        ] = f"{marker_names[group_index]}_pos"

    summary = {
        "project_dir": str(Path(project_dir).expanduser().resolve()),
        "output_dir": str(output_path),
        "threshold_source": threshold_source,
        "threshold_source_kind": threshold_source_kind,
        "generated_threshold_template": generated_threshold_template,
        "n_threshold_definitions": len(simple_classifiers),
        "n_marker_groups": len(matched_groups),
        "n_pos_columns": n_pos_columns,
        "pos_columns": active_pos_columns,
        "stale_celltype_columns": [
            column
            for column in replacement["removed_columns"]
            if column not in active_pos_columns
        ],
    }
    adata.uns["qxycell_thresholding"] = summary
    metadata["thresholding_applied"] = True
    metadata["threshold_source"] = threshold_source
    metadata["threshold_source_kind"] = threshold_source_kind
    metadata["generated_threshold_template"] = generated_threshold_template

    tables_dir = output_path / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    summary_path = tables_dir / "thresholding_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    stage_files = [summary_path]
    if generated_threshold_template:
        stage_files.append(Path(generated_threshold_template))
    complete_stage(
        adata,
        "thresholds",
        columns=active_pos_columns,
        files=stage_files,
        details={
            "source_kind": threshold_source_kind,
            "source_path": threshold_source,
            "n_pos_columns": n_pos_columns,
        },
    )
    if verbose:
        print(
            "Applied marker thresholds: "
            f"{n_pos_columns} positivity columns from {threshold_source}"
        )
    return summary


def threshold_from_classifiers(
    adata,
    project_dir: str | Path | None = None,
    *,
    output_dir: str | Path | None = None,
    image_col: str = "Image",
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply classifier JSON thresholds and save the applied threshold table.

    Supplied threshold tables are ignored. A successful run writes or replaces
    ``thresholds/classifier_thresholds.tsv`` in the active output directory.
    """

    metadata = getattr(adata, "uns", {}).get("qxycell", {})
    if project_dir is None and isinstance(metadata, dict):
        project_dir = metadata.get("project_dir")
    if project_dir is None:
        raise ValueError(
            "project_dir is required unless adata.uns['qxycell']['project_dir'] is set."
        )
    project_path = Path(project_dir).expanduser().resolve()
    classifiers = parse_classifiers(discover_classifier_files(project_path))
    conflicts = classifier_threshold_conflicts(classifiers)
    if conflicts:
        measurements = ", ".join(str(item["measurement_column"]) for item in conflicts)
        raise ValueError(
            "Conflicting object-classifier thresholds must be resolved in a threshold "
            "table before application: " + measurements
        )
    simple_paths = sorted(
        {str(item.path) for item in classifiers if item.is_simple}
    )
    if not simple_paths:
        raise ValueError("No usable QuPath classifier JSON thresholds are available.")
    output_path = _threshold_output_dir(adata, output_dir)
    return _apply_threshold_definitions(
        adata,
        classifiers,
        project_dir=project_path,
        output_path=output_path,
        threshold_source="|".join(simple_paths),
        threshold_source_kind="object_classifiers",
        threshold_table_writer=lambda: write_classifier_threshold_table(
            project_path,
            output_path,
            classifiers,
        ),
        image_col=image_col,
        verbose=verbose,
    )


def threshold_from_table(
    adata,
    threshold_file: str | Path,
    project_dir: str | Path | None = None,
    *,
    output_dir: str | Path | None = None,
    image_col: str = "Image",
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply only the named threshold table, ignoring classifier JSON thresholds."""

    metadata = getattr(adata, "uns", {}).get("qxycell", {})
    if project_dir is None and isinstance(metadata, dict):
        project_dir = metadata.get("project_dir")
    if project_dir is None:
        raise ValueError(
            "project_dir is required unless adata.uns['qxycell']['project_dir'] is set."
        )
    table_path = Path(threshold_file).expanduser().resolve()
    if not table_path.is_file():
        raise FileNotFoundError(f"Threshold table does not exist: {table_path}")
    unresolved = unresolved_threshold_conflicts(table_path)
    if unresolved:
        details = "; ".join(
            f"{item['marker'] or item['measurement_column']}: "
            + ", ".join(item["missing_images"])
            for item in unresolved
        )
        raise ValueError(
            "Threshold table contains unresolved classifier conflicts: " + details
        )
    classifiers = parse_threshold_files([table_path])
    return _apply_threshold_definitions(
        adata,
        classifiers,
        project_dir=project_dir,
        output_path=_threshold_output_dir(adata, output_dir),
        threshold_source=str(table_path),
        threshold_source_kind="manual_threshold_file",
        image_col=image_col,
        verbose=verbose,
    )


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


def import_cells(
    project_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    fail_on_check_error: bool = True,
    verbose: bool = True,
) -> Any:
    """Create the base AnnData cell checkpoint from QuPath measurement tables."""

    ad, np, pd = _import_runtime_dependencies()
    project_path = Path(project_dir).expanduser().resolve()
    if not project_path.exists():
        raise FileNotFoundError(f"QuPath project folder does not exist: {project_path}")
    output_path = resolve_output_dir(
        output_dir,
        project_dir=project_path,
        project_output_kind="run",
    )

    measurement_files = [
        summarize_measurement_file(path)
        for path in discover_measurement_files(project_path)
    ]
    messages = validate_measurement_files(measurement_files)
    errors = [message for message in messages if message.level == "error"]
    if fail_on_check_error and errors:
        raise RuntimeError("; ".join(message.message for message in errors))

    if verbose:
        print("Loading QuPath measurement table(s)...")
    measurements = _read_measurements(measurement_files, pd)
    for column in required_columns():
        if column not in measurements.columns:
            raise ValueError(f"Missing required QuPath measurement column: {column}")

    marker_columns = measurement_columns_for_threshold_template(measurement_files)
    if not marker_columns:
        raise ValueError("No mean/median measurement columns are available for adata.X import.")
    marker_names = _unique_marker_names_for_measurement_columns(marker_columns, [])
    missing_marker_columns = [
        column for column in marker_columns if column not in measurements.columns
    ]
    if missing_marker_columns:
        raise ValueError(
            "Measurement columns selected for adata.X are missing: "
            + ", ".join(str(column) for column in missing_marker_columns)
        )

    x = (
        measurements[marker_columns]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
        .to_numpy()
    )
    optional_obs_columns = ["TMA Core"] if "TMA Core" in measurements.columns else []
    obs_columns = (
        list(required_columns())
        + optional_obs_columns
        + ["quxy_source_file", "quxy_source_row"]
    )
    obs = measurements.loc[:, obs_columns].copy().rename(columns=CENTROID_OBS_COLUMN_RENAMES)
    obs["quxy_cell_id"] = obs["Image"].astype(str) + "::" + obs["Object ID"].astype(str)
    obs.index = obs["quxy_cell_id"].astype(str)
    spatial = obs[["Xµm", "Yµm"]].apply(pd.to_numeric, errors="coerce").to_numpy()
    if np.isnan(spatial).any():
        raise ValueError("Centroid X/Y columns contain missing or non-numeric values.")

    var = _build_var_dataframe_from_measurement_columns(marker_columns, marker_names, pd)
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.var_names = list(var.index)
    adata.obsm["spatial"] = spatial

    coreid_summary = None
    if "TMA Core" in adata.obs.columns:
        coreid_summary = assign_core_ids_from_measurements(adata, verbose=False)

    adata.uns["qxycell"] = {
        "project_dir": str(project_path),
        "output_dir": str(output_path),
        "run_dir": str(output_path),
        "tables_dir": str(output_path / "tables"),
        "n_measurement_files": len(measurement_files),
        "measurement_core_assignment": coreid_summary,
    }
    complete_stage(
        adata,
        "measurements",
        columns=list(adata.obs.columns),
        details={
            "project_dir": str(project_path),
            "measurement_files": [str(item.path) for item in measurement_files],
            "n_cells": int(adata.n_obs),
            "n_markers": int(adata.n_vars),
        },
    )
    if verbose:
        print(f"Created AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} markers")
    return adata


def add_annotations(
    adata,
    project_dir: str | Path | None = None,
    *,
    pixel_size_um: float = 0.28,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace GeoJSON-derived annotations and cell polygons on an AnnData object."""

    import pandas as pd

    pixel_size_um = _validate_pixel_size_um(pixel_size_um)
    metadata = getattr(adata, "uns", {}).get("qxycell", {})
    if project_dir is None and isinstance(metadata, dict):
        project_dir = metadata.get("project_dir")
    if project_dir is None:
        raise ValueError(
            "project_dir is required unless adata.uns['qxycell']['project_dir'] is set."
        )
    project_path = Path(project_dir).expanduser().resolve()
    geojson_files = summarize_geojson_files(discover_geojson_files(project_path))
    errors = [
        message
        for message in validate_geojson_files(geojson_files)
        if message.level == "error"
    ]
    if errors:
        raise RuntimeError("; ".join(message.message for message in errors))

    prepare_stage(adata, "annotations")
    adata.uns["qxycell_annotation_labels"] = {}
    columns_before = set(str(column) for column in adata.obs.columns)
    annotation_conflicts = _apply_annotations(
        adata,
        geojson_files,
        pixel_size_um=pixel_size_um,
        skip_annotation_labels=set(),
    )
    sample_annotation_cols = [
        column
        for column in adata.obs.columns
        if str(column).startswith("sample_annotation__")
    ]
    sample_summary = None
    if sample_annotation_cols:
        sample_summary = assign_samples(
            adata,
            annotation_prefix="sample_annotation__",
            sample_col="Sample",
            verbose=verbose,
        )
        adata.obs.drop(columns=sample_annotation_cols, inplace=True)
        label_map = adata.uns.get("qxycell_annotation_labels", {})
        if isinstance(label_map, dict):
            adata.uns["qxycell_annotation_labels"] = {
                column: label
                for column, label in label_map.items()
                if not str(column).startswith("sample_annotation__")
            }

    n_cell_polygons = _apply_cell_polygons(
        adata,
        geojson_files,
        pixel_size_um=pixel_size_um,
    )
    owned_columns = [
        str(column) for column in adata.obs.columns if str(column) not in columns_before
    ]
    label_map = adata.uns.get("qxycell_annotation_labels", {})
    assignment_rows = []
    if isinstance(label_map, dict):
        for column, label in label_map.items():
            if column not in adata.obs.columns:
                continue
            assignment_rows.append(
                {
                    "source_annotation": str(label),
                    "destination_column": str(column),
                    "n_assigned_cells": int(
                        adata.obs[column].fillna(False).astype(bool).sum()
                    ),
                }
            )

    output_path = resolve_output_dir(adata=adata)
    tables_dir = output_path / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    conflicts_path = tables_dir / "annotation_conflicts.csv"
    assignments_path = tables_dir / "annotation_assignments.csv"
    pd.DataFrame(annotation_conflicts).to_csv(conflicts_path, index=False)
    pd.DataFrame(assignment_rows).to_csv(assignments_path, index=False)
    summary = {
        "project_dir": str(project_path),
        "pixel_size_um": pixel_size_um,
        "n_geojson_files": len(geojson_files),
        "n_annotation_columns": sum(
            str(column).startswith("annotation__") for column in owned_columns
        ),
        "n_cell_polygons": int(n_cell_polygons),
        "sample_assignment": sample_summary,
    }
    adata.uns.setdefault("qxycell", {})["annotations"] = summary
    complete_stage(
        adata,
        "annotations",
        columns=owned_columns,
        files=[conflicts_path, assignments_path],
        details=summary,
    )
    if verbose:
        print(
            "Applied GeoJSON data: "
            f"{summary['n_annotation_columns']} annotation columns, "
            f"{n_cell_polygons:,} matched cell polygons"
        )
    return summary
