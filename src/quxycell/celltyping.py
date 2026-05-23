"""Cell type assignment from marker positivity columns."""

from __future__ import annotations

from pathlib import Path
from typing import Any


DEFAULT_CELLTYPE_DIR = Path("outputs") / "qxy_run" / "celltype"


def find_latest_celltype_yaml(celltype_dir: str | Path = DEFAULT_CELLTYPE_DIR) -> Path:
    """Return the newest YAML file in the cell type output folder."""

    celltype_dir = Path(celltype_dir).expanduser().resolve()
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


def apply_celltypes(
    adata,
    logic: str | Path | dict[str, Any] | None = None,
    *,
    celltype_column: str = "celltype",
    unknown_label: str = "Unknown",
    celltype_dir: str | Path = DEFAULT_CELLTYPE_DIR,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply ordered exclusive cell type rules and optional feature flags.

    Rules are evaluated in order. The first matching rule assigns the cell type.
    Functional feature flags are non-exclusive and are written to ``adata.obs``.
    """

    import numpy as np

    logic_source = None
    if logic is None:
        logic_source = find_latest_celltype_yaml(celltype_dir)
        logic = logic_source
    elif isinstance(logic, (str, Path)):
        logic_source = Path(logic).expanduser().resolve()

    if verbose and logic_source is not None:
        print(f"Using QUXYCell cell type logic YAML:\n{logic_source}")

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

    for rule in rules:
        name = str(rule.get("name") or "").strip()
        if not name:
            continue
        mask = _rule_mask(obs, rule, missing_references)
        assign_mask = mask & unassigned
        n_assigned = int(assign_mask.sum())
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
    }
    adata.uns["quxycell_celltyping"] = summary
    return summary
