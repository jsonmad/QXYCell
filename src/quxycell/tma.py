"""TMA core assignment from GeoJSON boundaries."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from quxycell.geojson import _classification_name
from quxycell.paths import resolve_output_dir
from quxycell.pipeline import _image_key
from quxycell.plotting import _resolve_spatial_key


def _safe_column(value: object) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(value).strip()).strip("_")
    return safe or "metadata"


def _geojson_paths(path: str | Path) -> list[Path]:
    path = Path(path).expanduser().resolve()
    if path.is_file():
        return [path]
    if path.is_dir():
        return sorted(
            [
                item
                for pattern in ("*.geojson", "*.json")
                for item in path.glob(pattern)
                if item.is_file()
            ]
        )
    raise FileNotFoundError(f"TMA GeoJSON path does not exist: {path}")


def _geojson_image_key(path: Path) -> str:
    stem = path.stem
    return stem[:-4] if stem.endswith(".ome") else stem


def _core_label(properties: dict[str, Any], fallback_index: int) -> str:
    for key in (
        "core_id",
        "Core ID",
        "core",
        "Core",
        "name",
        "Name",
        "label",
        "Label",
    ):
        value = properties.get(key)
        if value not in (None, ""):
            return str(value)
    label = _classification_name(properties)
    return label if label != "Unclassified" else f"core_{fallback_index:03d}"


def _scalar_metadata(properties: dict[str, Any]) -> dict[str, Any]:
    out = {}
    for key, value in properties.items():
        if isinstance(value, (str, int, float, bool)) or value is None:
            out[str(key)] = value
    return out


def _load_cores(core_geojson: str | Path, pixel_size_um: float):
    try:
        from shapely import affinity
        from shapely.geometry import shape
    except ImportError as exc:
        raise ImportError("TMA core assignment requires shapely.") from exc

    cores = []
    counter = 0
    for path in _geojson_paths(core_geojson):
        data = json.loads(path.read_text(errors="replace"))
        features = data.get("features", []) if isinstance(data, dict) else []
        image_key = _geojson_image_key(path)
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
            counter += 1
            cores.append(
                {
                    "core_label": _core_label(properties, counter),
                    "image_key": image_key,
                    "geometry": geometry,
                    "source": str(path),
                    "metadata": _scalar_metadata(properties),
                }
            )
    return cores


def _validate_no_core_overlaps(cores: list[dict[str, Any]]) -> None:
    overlap_rows = []
    for i, left in enumerate(cores):
        for right in cores[i + 1 :]:
            if left["image_key"] != right["image_key"]:
                continue
            if not left["geometry"].intersects(right["geometry"]):
                continue
            intersection = left["geometry"].intersection(right["geometry"])
            if intersection.area > 0:
                overlap_rows.append(
                    f"{left['core_label']} overlaps {right['core_label']} "
                    f"in {left['image_key']} (area={intersection.area:.3f})"
                )
    if overlap_rows:
        raise ValueError(
            "Overlapping TMA cores are not allowed:\n" + "\n".join(overlap_rows[:20])
        )


def _cores_for_sample(
    cores: list[dict[str, Any]],
    sample: str,
    *,
    n_samples: int,
) -> list[dict[str, Any]]:
    sample_key = _image_key(sample)
    exact = [core for core in cores if core["image_key"] == sample_key]
    if exact:
        return exact
    image_keys = sorted({core["image_key"] for core in cores})
    if len(image_keys) == 1 and n_samples == 1:
        return cores
    return []


def assign_tma_cores(
    adata,
    core_geojson: str | Path,
    *,
    sample_col: str = "Image",
    spatial_key: str | None = None,
    pixel_size_um: float = 1.0,
    core_col: str = "tma_core",
    metadata_prefix: str = "tma_",
    unassigned_label: str = "Unassigned",
    output_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign cells to non-overlapping TMA cores from GeoJSON boundaries.

    TMA GeoJSON files are matched to ``adata.obs[sample_col]`` using the
    GeoJSON filename stem. With the default ``sample_col="Image"``, this matches
    QuPath image names. A single unmatched GeoJSON is only applied globally when
    the AnnData object contains one sample.
    """

    import pandas as pd
    from shapely.geometry import Point

    if sample_col not in adata.obs.columns:
        raise KeyError(f"sample_col not found in adata.obs: {sample_col}")

    spatial_key = _resolve_spatial_key(adata, spatial_key)
    cores = _load_cores(core_geojson, pixel_size_um=pixel_size_um)
    if not cores:
        raise ValueError(f"No polygon TMA cores found in {core_geojson}")
    _validate_no_core_overlaps(cores)

    metadata_keys = sorted({key for core in cores for key in core["metadata"]})
    metadata_columns = {key: metadata_prefix + _safe_column(key) for key in metadata_keys}

    adata.obs[core_col] = unassigned_label
    adata.obs[f"{core_col}_source"] = ""
    for column in metadata_columns.values():
        adata.obs[column] = pd.NA

    coords = adata.obsm[spatial_key]
    sample_values = adata.obs[sample_col].astype(str)
    n_samples = int(sample_values.nunique())
    assigned = 0
    unmatched_samples = set()

    for sample in sorted(sample_values.unique()):
        sample_cores = _cores_for_sample(cores, sample, n_samples=n_samples)
        if not sample_cores:
            unmatched_samples.add(sample)
            continue
        sample_mask = sample_values == sample
        positions = sample_mask.to_numpy().nonzero()[0]
        for position in positions:
            point = Point(float(coords[position, 0]), float(coords[position, 1]))
            hits = [core for core in sample_cores if core["geometry"].covers(point)]
            if len(hits) > 1:
                labels = ", ".join(core["core_label"] for core in hits)
                raise ValueError(
                    "Cell assigned to multiple TMA cores. "
                    f"sample={sample}, obs_index={adata.obs.index[position]}, cores={labels}"
                )
            if not hits:
                continue
            core = hits[0]
            obs_index = adata.obs.index[position]
            adata.obs.at[obs_index, core_col] = core["core_label"]
            adata.obs.at[obs_index, f"{core_col}_source"] = core["source"]
            for key, column in metadata_columns.items():
                adata.obs.at[obs_index, column] = core["metadata"].get(key, pd.NA)
            assigned += 1

    out_dir = resolve_output_dir(output_dir, adata=adata) / "tma"
    out_dir.mkdir(parents=True, exist_ok=True)
    core_table = pd.DataFrame(
        [
            {
                "core_label": core["core_label"],
                "image_key": core["image_key"],
                "source": core["source"],
                "area": float(core["geometry"].area),
                **core["metadata"],
            }
            for core in cores
        ]
    )
    core_table_path = out_dir / "tma_cores.tsv"
    core_table.to_csv(core_table_path, sep="\t", index=False)

    counts = (
        adata.obs[core_col].astype(str)
        .value_counts(dropna=False)
        .rename_axis(core_col)
        .reset_index(name="n_cells")
    )
    counts_path = out_dir / "tma_core_counts.tsv"
    counts.to_csv(counts_path, sep="\t", index=False)

    summary = {
        "core_col": core_col,
        "n_cores": len(cores),
        "n_cells": int(adata.n_obs),
        "n_assigned_cells": int(assigned),
        "n_unassigned_cells": int(adata.n_obs - assigned),
        "sample_col": sample_col,
        "core_image_keys": sorted({core["image_key"] for core in cores}),
        "unmatched_samples": sorted(unmatched_samples),
        "core_table_tsv": str(core_table_path),
        "counts_tsv": str(counts_path),
        "metadata_columns": list(metadata_columns.values()),
    }
    adata.uns["quxycell_tma"] = summary

    if verbose:
        print(
            "Assigned TMA cores: "
            f"{assigned:,}/{adata.n_obs:,} cells across {len(cores)} cores"
        )
        print(f"Saved TMA summary:\n{counts_path}")

    return summary
