"""Sample metadata import helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from qxycell.paths import resolve_output_dir


def _read_metadata_table(metadata: str | Path | Any):
    import pandas as pd

    if hasattr(metadata, "copy") and hasattr(metadata, "columns"):
        return metadata.copy()

    path = Path(metadata).expanduser().resolve()
    if path.suffix.lower() in {".tsv", ".txt"}:
        sep = "\t"
    else:
        sep = ","
    return pd.read_csv(path, sep=sep)


def _h5ad_safe_metadata_series(values, source_column):
    import pandas as pd
    from pandas.api.types import is_bool_dtype, is_numeric_dtype

    if is_bool_dtype(source_column):
        return values.astype("boolean")
    if is_numeric_dtype(source_column):
        return pd.to_numeric(values, errors="coerce")
    return values.astype("string").astype("category")


def _metadata_key_series(values):
    keys = values.astype("string").str.strip()
    return keys.fillna("nan")


def add_metadata(
    adata,
    metadata: str | Path | Any,
    *,
    sample_col: str = "Image",
    metadata_sample_col: str | None = None,
    columns: list[str] | None = None,
    prefix: str = "",
    overwrite: bool = False,
    output_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add sample-level metadata columns to ``adata.obs``.

    ``metadata`` can be a CSV/TSV path or a pandas DataFrame. One row per sample
    is required. ``sample_col`` can be any column in ``adata.obs``; the default
    is the QuPath ``Image`` column.
    """

    import pandas as pd

    if sample_col not in adata.obs.columns:
        raise KeyError(f"sample_col not found in adata.obs: {sample_col}")

    table = _read_metadata_table(metadata)
    if metadata_sample_col is None:
        metadata_sample_col = sample_col if sample_col in table.columns else str(table.columns[0])
    if metadata_sample_col not in table.columns:
        raise KeyError(f"metadata_sample_col not found in metadata: {metadata_sample_col}")

    table = table.copy()
    table[metadata_sample_col] = _metadata_key_series(table[metadata_sample_col])
    duplicated = table[table[metadata_sample_col].duplicated()][metadata_sample_col].unique().tolist()
    if duplicated:
        raise ValueError(f"Sample metadata contains duplicate keys: {duplicated}")

    if columns is None:
        columns = [str(column) for column in table.columns if column != metadata_sample_col]
    missing_columns = [column for column in columns if column not in table.columns]
    if missing_columns:
        raise KeyError(f"Metadata columns not found: {missing_columns}")

    obs_samples = _metadata_key_series(adata.obs[sample_col])
    metadata_indexed = table.set_index(metadata_sample_col, drop=False)
    metadata_samples = set(metadata_indexed.index.astype(str))
    adata_samples = set(obs_samples.unique())
    matched_samples = sorted(adata_samples & metadata_samples)
    missing_in_metadata = sorted(adata_samples - metadata_samples)
    unused_metadata = sorted(metadata_samples - adata_samples)

    added_columns = []
    for column in columns:
        output_column = f"{prefix}{column}"
        if output_column in adata.obs.columns and not overwrite:
            raise ValueError(
                f"Column already exists in adata.obs: {output_column}. "
                "Use overwrite=True to replace it."
            )
        mapping = metadata_indexed[column].to_dict()
        values = obs_samples.map(mapping)
        adata.obs[output_column] = _h5ad_safe_metadata_series(values, table[column])
        added_columns.append(output_column)

    summary = {
        "sample_col": sample_col,
        "metadata_sample_col": metadata_sample_col,
        "n_adata_samples": len(adata_samples),
        "n_metadata_samples": len(metadata_samples),
        "n_matched_samples": len(matched_samples),
        "missing_in_metadata": missing_in_metadata,
        "unused_metadata": unused_metadata,
        "added_columns": added_columns,
    }
    adata.uns["qxycell_sample_metadata"] = summary

    out_dir = resolve_output_dir(output_dir, adata=adata) / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "sample_metadata_summary.tsv"
    applied_path = out_dir / "sample_metadata_applied.tsv"
    pd.DataFrame(
        [
            {"metric": key, "value": value if not isinstance(value, list) else "; ".join(map(str, value))}
            for key, value in summary.items()
        ]
    ).to_csv(summary_path, sep="\t", index=False)
    table.to_csv(applied_path, sep="\t", index=False)
    summary["summary_tsv"] = str(summary_path)
    summary["metadata_tsv"] = str(applied_path)

    if verbose:
        print(
            "Added sample metadata: "
            f"{len(added_columns)} columns, {len(matched_samples)}/{len(adata_samples)} samples matched"
        )
        print(f"Saved sample metadata summary:\n{summary_path}")

    return summary
