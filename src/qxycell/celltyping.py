"""Cell type assignment from marker positivity columns."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from qxycell.paths import latest_timestamped_output_dir
from qxycell.paths import output_dir_from_adata
from qxycell.paths import resolve_output_dir


def _safe_name(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return safe.strip("_") or "value"


def _default_celltype_dir(adata=None, logic_source: Path | None = None) -> Path:
    adata_output_dir = output_dir_from_adata(adata) if adata is not None else None
    if adata_output_dir is not None:
        return adata_output_dir / "celltype"
    if logic_source is not None:
        return logic_source.parent
    latest_output_dir = latest_timestamped_output_dir()
    if latest_output_dir is not None:
        return latest_output_dir / "celltype"
    return resolve_output_dir() / "celltype"


def find_latest_celltype_yaml(celltype_dir: str | Path | None = None) -> Path:
    """Return the newest YAML file in the cell type output folder."""

    celltype_dir = (
        Path(celltype_dir).expanduser().resolve()
        if celltype_dir is not None
        else _default_celltype_dir()
    )
    candidates = [
        path
        for pattern in ("*.yaml", "*.yml")
        for path in celltype_dir.glob(pattern)
        if path.is_file()
    ]
    if not candidates:
        raise FileNotFoundError(
            "No cell type YAML files were found in "
            f"{celltype_dir}. Provide a YAML path, or save one there first."
        )
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _load_yaml(path: str | Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("Cell type logic YAML support requires PyYAML.") from exc

    path = Path(path).expanduser().resolve()
    if path.suffix.lower() not in {".yaml", ".yml"}:
        raise ValueError(
            "Cell type logic must be a .yaml or .yml file. "
            f"Got: {path}"
        )
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Cell type logic must be a YAML mapping: {path}")
    return data


def load_celltype_logic(logic: str | Path | dict[str, Any]) -> dict[str, Any]:
    """Load cell type logic from a YAML file or mapping.

    Expected shape:

    ```yaml
    rules:
      - name: CD8 T cells
        positive: [CD3, CD8]
        negative: [CD4]
    features:
      PD1_positive:
        positive: [PD1]
    derived_features:
      PD1_GZMB_double_positive:
        all_of: [PD1_positive, GZMB_positive]
    ```
    """

    if isinstance(logic, (str, Path)):
        logic = _load_yaml(logic)
    if not isinstance(logic, dict):
        raise TypeError("celltype logic must be a YAML path or dictionary.")

    rules = logic.get("rules", [])
    features = logic.get("features", {})
    derived_features = logic.get("derived_features", {})
    if rules is None:
        rules = []
    if features is None:
        features = {}
    if derived_features is None:
        derived_features = {}
    if not isinstance(rules, list):
        raise ValueError("celltype logic 'rules' must be a list.")
    if not isinstance(features, dict):
        raise ValueError("celltype logic 'features' must be a mapping.")
    if not isinstance(derived_features, dict):
        raise ValueError("celltype logic 'derived_features' must be a mapping.")

    return {
        "rules": rules,
        "features": features,
        "derived_features": derived_features,
    }


def _obs_bool_array(obs, column_name: str):
    import numpy as np

    if column_name not in obs.columns:
        return np.zeros(len(obs), dtype=bool)
    return obs[column_name].fillna(0).astype(bool).to_numpy()


def _resolve_reference_column(obs, reference_name: str) -> str | None:
    reference_name = str(reference_name)
    if reference_name in obs.columns:
        return reference_name
    pos_column = f"{reference_name}_pos"
    if pos_column in obs.columns:
        return pos_column
    return None


def _rule_mask(obs, rule: dict[str, Any], missing_references: set[str]):
    import numpy as np

    mask = np.ones(len(obs), dtype=bool)

    for reference_name in rule.get("positive", []) or []:
        column = _resolve_reference_column(obs, str(reference_name))
        if column is None:
            missing_references.add(str(reference_name))
            return np.zeros(len(obs), dtype=bool)
        mask &= _obs_bool_array(obs, column)

    for reference_name in rule.get("negative", []) or []:
        column = _resolve_reference_column(obs, str(reference_name))
        if column is None:
            continue
        mask &= ~_obs_bool_array(obs, column)

    any_positive = rule.get("any_positive", []) or []
    if any_positive:
        any_mask = np.zeros(len(obs), dtype=bool)
        any_found = False
        for reference_name in any_positive:
            column = _resolve_reference_column(obs, str(reference_name))
            if column is None:
                missing_references.add(str(reference_name))
                continue
            any_found = True
            any_mask |= _obs_bool_array(obs, column)
        mask &= any_mask if any_found else np.zeros(len(obs), dtype=bool)

    return mask


def _derived_mask(obs, spec: dict[str, Any], missing_references: set[str]):
    import numpy as np

    mask = np.ones(len(obs), dtype=bool)

    for key in ("positive", "all_of"):
        for reference_name in spec.get(key, []) or []:
            column = _resolve_reference_column(obs, str(reference_name))
            if column is None:
                missing_references.add(str(reference_name))
                return np.zeros(len(obs), dtype=bool)
            mask &= _obs_bool_array(obs, column)

    for key in ("negative", "none_of"):
        for reference_name in spec.get(key, []) or []:
            column = _resolve_reference_column(obs, str(reference_name))
            if column is None:
                continue
            mask &= ~_obs_bool_array(obs, column)

    any_of = spec.get("any_of", []) or []
    if any_of:
        any_mask = np.zeros(len(obs), dtype=bool)
        any_found = False
        for reference_name in any_of:
            column = _resolve_reference_column(obs, str(reference_name))
            if column is None:
                missing_references.add(str(reference_name))
                continue
            any_found = True
            any_mask |= _obs_bool_array(obs, column)
        mask &= any_mask if any_found else np.zeros(len(obs), dtype=bool)

    return mask


def _stringify_definition_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple, set)):
        return ", ".join(str(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True)
    return str(value)


def _rule_definition(rule: dict[str, Any]) -> str:
    parts = []
    for key, label in (
        ("positive", "positive"),
        ("negative", "negative"),
        ("any_positive", "any_positive"),
    ):
        value = rule.get(key)
        if value:
            parts.append(f"{label}: {_stringify_definition_value(value)}")
    return "; ".join(parts)


def _write_rule_summary_table(
    rule_diagnostics: list[dict[str, Any]],
    *,
    celltype_dir: str | Path,
    logic_source: Path | None,
    celltype_column: str,
    unknown_label: str,
    n_cells: int,
) -> Path:
    import pandas as pd

    celltype_dir = Path(celltype_dir).expanduser().resolve()
    celltype_dir.mkdir(parents=True, exist_ok=True)
    logic_name = _safe_name(logic_source.stem) if logic_source is not None else "dict"
    output_path = celltype_dir / f"celltype_rules_summary_{logic_name}.tsv"

    known_keys = {"name", "positive", "negative", "any_positive"}
    rows = []
    for diagnostic in rule_diagnostics:
        rule = diagnostic["rule"]
        index = diagnostic["rule_order"]
        name = str(rule.get("name") or "").strip()
        extra = {key: value for key, value in rule.items() if key not in known_keys}
        rows.append(
            {
                "rule_order": index,
                "celltype": name,
                "positive": _stringify_definition_value(rule.get("positive")),
                "negative": _stringify_definition_value(rule.get("negative")),
                "any_positive": _stringify_definition_value(rule.get("any_positive")),
                "definition": _rule_definition(rule),
                "raw_matching_cells": int(diagnostic["raw_matching_cells"]),
                "assigned_cells": int(diagnostic["assigned_cells"]),
                "blocked_by_prior_rules": int(diagnostic["blocked_by_prior_rules"]),
                "overlap_with_other_rules_cells": int(diagnostic["overlap_with_other_rules_cells"]),
                "raw_fraction": int(diagnostic["raw_matching_cells"]) / max(n_cells, 1),
                "assigned_fraction": int(diagnostic["assigned_cells"]) / max(n_cells, 1),
                "missing_references": ", ".join(diagnostic["missing_references"]),
                "celltype_column": celltype_column,
                "extra_definition": json.dumps(extra, sort_keys=True) if extra else "",
            }
        )

    table = pd.DataFrame(
        rows,
        columns=[
            "rule_order",
            "celltype",
            "positive",
            "negative",
            "any_positive",
            "definition",
            "raw_matching_cells",
            "assigned_cells",
            "blocked_by_prior_rules",
            "overlap_with_other_rules_cells",
            "raw_fraction",
            "assigned_fraction",
            "missing_references",
            "celltype_column",
            "extra_definition",
        ],
    )
    table.to_csv(output_path, sep="\t", index=False)
    return output_path


def apply_celltypes(
    adata,
    logic: str | Path | dict[str, Any] | None = None,
    *,
    celltype_column: str = "celltype",
    unknown_label: str = "Unknown",
    celltype_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply ordered exclusive cell type rules and optional feature flags.

    Rules are evaluated in order. The first matching rule assigns the cell type.
    Functional feature flags are non-exclusive and are written to ``adata.obs``.
    """

    import numpy as np

    logic_source = None
    if isinstance(logic, (str, Path)):
        logic_source = Path(logic).expanduser().resolve()
    resolved_celltype_dir = (
        Path(celltype_dir).expanduser().resolve()
        if celltype_dir is not None
        else _default_celltype_dir(adata, logic_source)
    )
    if logic is None:
        logic_source = find_latest_celltype_yaml(resolved_celltype_dir)
        logic = logic_source

    if verbose and logic_source is not None:
        print(f"Using QXYCell cell type logic YAML:\n{logic_source}")

    bundle = load_celltype_logic(logic)
    obs = adata.obs
    rules = [rule for rule in bundle["rules"] if isinstance(rule, dict)]
    features = {
        str(name): spec
        for name, spec in bundle["features"].items()
        if isinstance(spec, dict)
    }
    derived_features = {
        str(name): spec
        for name, spec in bundle["derived_features"].items()
        if isinstance(spec, dict)
    }

    celltypes = np.full(len(obs), unknown_label, dtype=object)
    unassigned = np.ones(len(obs), dtype=bool)
    missing_references: set[str] = set()
    assigned_counts: dict[str, int] = {}
    rule_diagnostics: list[dict[str, Any]] = []
    rule_masks = []

    for rule_order, rule in enumerate(rules, start=1):
        name = str(rule.get("name") or "").strip()
        if not name:
            continue
        rule_missing: set[str] = set()
        mask = _rule_mask(obs, rule, rule_missing)
        missing_references.update(rule_missing)
        rule_masks.append(mask)
        rule_diagnostics.append(
            {
                "rule_order": rule_order,
                "rule": rule,
                "name": name,
                "mask": mask,
                "raw_matching_cells": int(mask.sum()),
                "assigned_cells": 0,
                "blocked_by_prior_rules": 0,
                "overlap_with_other_rules_cells": 0,
                "missing_references": sorted(rule_missing),
            }
        )

    if rule_masks:
        match_counts = np.vstack(rule_masks).sum(axis=0)
        multi_rule_match_count = int((match_counts > 1).sum())
    else:
        match_counts = np.zeros(len(obs), dtype=int)
        multi_rule_match_count = 0

    for diagnostic in rule_diagnostics:
        name = diagnostic["name"]
        mask = diagnostic["mask"]
        assign_mask = mask & unassigned
        n_assigned = int(assign_mask.sum())
        diagnostic["assigned_cells"] = n_assigned
        diagnostic["blocked_by_prior_rules"] = int(mask.sum()) - n_assigned
        diagnostic["overlap_with_other_rules_cells"] = int((mask & (match_counts > 1)).sum())
        if n_assigned:
            celltypes[assign_mask] = name
            unassigned[assign_mask] = False
            assigned_counts[name] = n_assigned

    obs[celltype_column] = celltypes

    feature_columns = []
    derived_feature_columns = []
    for feature_name, spec in features.items():
        column = str(feature_name)
        obs[column] = _rule_mask(obs, spec, missing_references).astype("int8")
        feature_columns.append(column)

    for feature_name, spec in derived_features.items():
        column = str(feature_name)
        obs[column] = _derived_mask(obs, spec, missing_references).astype("int8")
        derived_feature_columns.append(column)

    summary = {
        "celltype_column": celltype_column,
        "unknown_label": unknown_label,
        "logic_source": str(logic_source) if logic_source is not None else "<dict>",
        "n_rules": len(rules),
        "n_features": len(feature_columns),
        "n_derived_features": len(derived_feature_columns),
        "feature_columns": feature_columns,
        "derived_feature_columns": derived_feature_columns,
        "missing_references": sorted(missing_references),
        "assigned_counts": assigned_counts,
        "unknown_count": int((obs[celltype_column] == unknown_label).sum()),
        "zero_raw_match_rules": [
            item["name"] for item in rule_diagnostics if item["raw_matching_cells"] == 0
        ],
        "zero_assigned_rules": [
            item["name"] for item in rule_diagnostics if item["assigned_cells"] == 0
        ],
        "multi_rule_match_count": multi_rule_match_count,
    }
    rule_summary_tsv = _write_rule_summary_table(
        rule_diagnostics,
        celltype_dir=resolved_celltype_dir,
        logic_source=logic_source,
        celltype_column=celltype_column,
        unknown_label=unknown_label,
        n_cells=len(obs),
    )
    summary["rule_summary_tsv"] = str(rule_summary_tsv)
    adata.uns["qxycell_celltyping"] = summary
    if verbose:
        print(f"Saved cell type rule summary TSV:\n{rule_summary_tsv}")
    return summary
