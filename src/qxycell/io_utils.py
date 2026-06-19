"""Save and load helpers for QXYCell AnnData outputs."""

from __future__ import annotations

from pathlib import Path

from qxycell.paths import _is_legacy_default_output_path
from qxycell.paths import latest_timestamped_output_dir
from qxycell.paths import resolve_output_dir


def _has_legacy_parent(path: Path) -> bool:
    return any(_is_legacy_default_output_path(parent) for parent in path.parents)


def _h5ad_filename(output_path: Path) -> str:
    """Return the h5ad filename for a given output directory, including timestamp."""
    folder = output_path.name
    ts = folder.removeprefix("qxy_outputs_") if folder.startswith("qxy_outputs_") else ""
    return f"qxycell_{ts}.h5ad" if ts else "qxycell.h5ad"


def _default_h5ad_path(adata=None, output_dir: str | Path | None = None) -> Path:
    if output_dir is None and adata is not None:
        metadata = getattr(adata, "uns", {}).get("qxycell", {})
        if isinstance(metadata, dict) and metadata.get("h5ad_path"):
            h5ad_path = Path(metadata["h5ad_path"]).expanduser().resolve()
            if not _has_legacy_parent(h5ad_path):
                return h5ad_path
    output_path = resolve_output_dir(output_dir, adata=adata)
    return output_path / "run" / "h5ad" / _h5ad_filename(output_path)


def _prepare_obs_for_h5ad(adata) -> None:
    index_name = getattr(adata.obs.index, "name", None)
    if index_name is None or index_name not in adata.obs.columns:
        return

    index_values = adata.obs.index.astype(str)
    column_values = adata.obs[index_name].astype(str).to_numpy()
    if (index_values == column_values).all():
        adata.obs.index.name = None
        return

    raise ValueError(
        f"adata.obs.index.name {index_name!r} also exists as an obs column, "
        "but the index values do not match that column. Rename the index or "
        "column before saving."
    )


def save(
    adata,
    path: str | Path | None = None,
    *,
    output_dir: str | Path | None = None,
    verbose: bool = True,
) -> Path:
    """Save an AnnData object to the current QXYCell H5AD path."""

    h5ad_path = Path(path).expanduser().resolve() if path is not None else _default_h5ad_path(
        adata,
        output_dir,
    )
    h5ad_path.parent.mkdir(parents=True, exist_ok=True)
    _prepare_obs_for_h5ad(adata)
    adata.write_h5ad(h5ad_path)

    output_path = h5ad_path.parents[2] if h5ad_path.match("*/run/h5ad/*.h5ad") else h5ad_path.parent
    metadata = adata.uns.setdefault("qxycell", {})
    metadata["output_dir"] = str(output_path)
    metadata["run_dir"] = str(output_path / "run")
    metadata["h5ad_path"] = str(h5ad_path)
    metadata["tables_dir"] = str(output_path / "run" / "tables")

    if verbose:
        print(f"Saved QXYCell H5AD:\n{h5ad_path}")
    return h5ad_path


def _resolve_h5ad_input(path_or_output_dir: str | Path) -> Path:
    path = Path(path_or_output_dir).expanduser().resolve()
    if path.is_file():
        return path
    h5ad_dir = path / "run" / "h5ad"
    # Try timestamped name first, then fall back to any qxycell*.h5ad in the folder
    candidates = sorted(h5ad_dir.glob("qxycell*.h5ad")) if h5ad_dir.is_dir() else []
    if candidates:
        return max(candidates, key=lambda p: p.stat().st_mtime)
    raise FileNotFoundError(f"No QXYCell H5AD found in {h5ad_dir}")


def load(path_or_output_dir: str | Path):
    """Load a QXYCell H5AD from a file path or output directory."""

    import anndata as ad

    h5ad_path = _resolve_h5ad_input(path_or_output_dir)
    return ad.read_h5ad(h5ad_path)


def load_latest(base_dir: str | Path = "."):
    """Load the latest timestamped QXYCell output in ``base_dir``."""

    latest = latest_timestamped_output_dir(base_dir)
    if latest is None:
        raise FileNotFoundError(f"No qxy_outputs_YYMMDD-HHMM folders found in {Path(base_dir).resolve()}")
    return load(latest)
