import pytest

from qxycell.celltyping import load_celltype_logic


def test_celltype_logic_rejects_notebook_path(tmp_path):
    notebook_path = tmp_path / "celltype_logic.ipynb"
    notebook_path.write_text('{"cells": [], "metadata": {}}', encoding="utf-8")

    with pytest.raises(ValueError, match=r"\.yaml or \.yml"):
        load_celltype_logic(notebook_path)


def test_celltype_logic_loads_yaml_rules(tmp_path):
    yaml_path = tmp_path / "celltype_logic.yaml"
    yaml_path.write_text(
        "rules:\n"
        "  - name: T cell\n"
        "    positive: [CD3]\n",
        encoding="utf-8",
    )

    logic = load_celltype_logic(yaml_path)

    assert logic["rules"][0]["name"] == "T cell"
