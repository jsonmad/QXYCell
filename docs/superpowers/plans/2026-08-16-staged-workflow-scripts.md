# Staged Workflow Scripts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add public, independently runnable Python examples for every QXYCell staged-workflow step.

**Architecture:** A shared `config.py` supplies paths and settings to six small scripts. Stage 1 creates one fixed output checkpoint; every later script reloads that checkpoint and invokes exactly one staged API function. A local ignored structural test verifies the public examples without publishing tests.

**Tech Stack:** Python 3.10+, `pathlib`, QXYCell public API, pytest, Ruff, static HTML.

## Global Constraints

- Store public examples under `examples/staged_workflow/`.
- Include no private paths, datasets, email addresses, or research-specific values.
- Keep classifier-only and table-only thresholding in separate scripts with no fallback.
- Let QXYCell stage functions perform checkpointing; do not add duplicate save logic.
- Keep tests local and ignored; do not publish them.
- Keep examples outside the wheel and source distribution.

---

### Task 1: Shared configuration and staged scripts

**Files:**
- Create: `tests/test_staged_workflow_examples.py` (local and ignored)
- Create: `examples/staged_workflow/AGENTS.md`
- Create: `examples/staged_workflow/config.py`
- Create: `examples/staged_workflow/01_import_measurements.py`
- Create: `examples/staged_workflow/02_add_annotations.py`
- Create: `examples/staged_workflow/03a_threshold_from_classifiers.py`
- Create: `examples/staged_workflow/03b_threshold_from_table.py`
- Create: `examples/staged_workflow/04_generate_celltype_prompt.py`
- Create: `examples/staged_workflow/05_apply_celltypes.py`

**Interfaces:**
- Consumes: `qxy.import_measurements`, `qxy.load`, `qxy.add_annotations`, `qxy.threshold_from_classifiers`, `qxy.threshold_from_table`, `qxy.celltype_prompt`, and `qxy.celltype`.
- Produces: six import-safe scripts sharing `PROJECT_DIR`, `OUTPUT_DIR`, `THRESHOLD_TABLE`, `CELLTYPE_YAML`, `PIXEL_SIZE_UM`, and `CELLTYPE_CONTEXT` from `config.py`.

- [ ] **Step 1: Write the failing structural test**

Create a local ignored pytest that parses each script with `ast`, verifies its expected QXYCell calls, verifies stages 2–5 load `OUTPUT_DIR`, and verifies the threshold scripts do not call the other threshold route.

```python
from __future__ import annotations

import ast
from pathlib import Path


EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "staged_workflow"
EXPECTED_CALLS = {
    "01_import_measurements.py": {"import_measurements"},
    "02_add_annotations.py": {"load", "add_annotations"},
    "03a_threshold_from_classifiers.py": {"load", "threshold_from_classifiers"},
    "03b_threshold_from_table.py": {"load", "threshold_from_table"},
    "04_generate_celltype_prompt.py": {"load", "celltype_prompt"},
    "05_apply_celltypes.py": {"load", "celltype"},
}


def _qxy_calls(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return {
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "qxy"
    }


def test_staged_examples_have_one_explicit_stage_each():
    for filename, expected in EXPECTED_CALLS.items():
        assert _qxy_calls(EXAMPLE_DIR / filename) == expected


def test_threshold_examples_do_not_fall_back():
    classifier_calls = _qxy_calls(EXAMPLE_DIR / "03a_threshold_from_classifiers.py")
    table_calls = _qxy_calls(EXAMPLE_DIR / "03b_threshold_from_table.py")
    assert "threshold_from_table" not in classifier_calls
    assert "threshold_from_classifiers" not in table_calls


def test_config_exposes_shared_settings():
    source = (EXAMPLE_DIR / "config.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    names = {
        target.id
        for node in tree.body
        if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Name)
    }
    assert {
        "PROJECT_DIR",
        "OUTPUT_DIR",
        "THRESHOLD_TABLE",
        "CELLTYPE_YAML",
        "PIXEL_SIZE_UM",
        "CELLTYPE_CONTEXT",
    } <= names
```

- [ ] **Step 2: Run the structural test and verify RED**

Run:

```bash
python -m pytest tests/test_staged_workflow_examples.py -q
```

Expected: FAIL because `examples/staged_workflow/*.py` do not exist.

- [ ] **Step 3: Add the local subtree contract and shared configuration**

`AGENTS.md` will require one stage per script, shared configuration, no private paths, no implicit threshold selection, and import-safe `main()` guards. `config.py` will contain:

```python
"""Edit these settings once before running the staged workflow examples."""

from pathlib import Path


PROJECT_DIR = Path("/path/to/qupath_project")
OUTPUT_DIR = Path("/path/to/qxycell_output")
THRESHOLD_TABLE = PROJECT_DIR / "thresholds.tsv"
CELLTYPE_YAML = OUTPUT_DIR / "celltype" / "celltype_logic.yaml"
PIXEL_SIZE_UM = 0.28
CELLTYPE_CONTEXT = (
    "Describe the tissue, disease, experimental groups, and expected cell "
    "populations here. The generated rules require expert review."
)
```

- [ ] **Step 4: Implement the six minimal scripts**

Use this pattern for stage 1:

```python
"""Stage 1: create the base AnnData checkpoint from QuPath measurements."""

import qxycell as qxy

from config import OUTPUT_DIR, PROJECT_DIR


def main() -> None:
    qxy.import_measurements(PROJECT_DIR, output_dir=OUTPUT_DIR)


if __name__ == "__main__":
    main()
```

Implement stages 2 through 5 as follows:

```python
# 02_add_annotations.py
import qxycell as qxy
from config import OUTPUT_DIR, PIXEL_SIZE_UM, PROJECT_DIR

def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.add_annotations(
        adata, project_dir=PROJECT_DIR, pixel_size_um=PIXEL_SIZE_UM
    )

if __name__ == "__main__":
    main()
```

```python
# 03a_threshold_from_classifiers.py
import qxycell as qxy
from config import OUTPUT_DIR, PROJECT_DIR

def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.threshold_from_classifiers(adata, project_dir=PROJECT_DIR)

if __name__ == "__main__":
    main()
```

```python
# 03b_threshold_from_table.py
import qxycell as qxy
from config import OUTPUT_DIR, PROJECT_DIR, THRESHOLD_TABLE

def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.threshold_from_table(
        adata, THRESHOLD_TABLE, project_dir=PROJECT_DIR
    )

if __name__ == "__main__":
    main()
```

```python
# 04_generate_celltype_prompt.py
import qxycell as qxy
from config import CELLTYPE_CONTEXT, OUTPUT_DIR

def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.celltype_prompt(adata, context=CELLTYPE_CONTEXT)

if __name__ == "__main__":
    main()
```

```python
# 05_apply_celltypes.py
import qxycell as qxy
from config import CELLTYPE_YAML, OUTPUT_DIR

def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.celltype(adata, CELLTYPE_YAML)

if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run the structural test and verify GREEN**

Run:

```bash
python -m pytest tests/test_staged_workflow_examples.py -q
```

Expected: `3 passed`.

- [ ] **Step 6: Compile, lint, and import the scripts without execution**

Run:

```bash
python -m compileall -q examples/staged_workflow
ruff check examples/staged_workflow
```

Then load each numbered script with `runpy.run_path(..., run_name="example_import_check")`. Expected: no workflow stage executes and no exception is raised.

- [ ] **Step 7: Commit the public scripts**

Stage only `examples/staged_workflow/`; verify the ignored test is absent from the index; commit with:

```bash
git commit -m "docs: add staged workflow scripts"
```

---

### Task 2: Example guide and documentation links

**Files:**
- Create: `examples/staged_workflow/README.md`
- Modify: `README.md`
- Modify: `docs/index.html`
- Modify: `docs/QXYCell_overview.html`

**Interfaces:**
- Consumes: the seven public example files from Task 1.
- Produces: a beginner-facing execution guide and discoverable links from the main documentation.

- [ ] **Step 1: Extend the local structural test and verify RED**

Add assertions that `README.md` exists in the example directory, names all six scripts, says to choose exactly one stage 3 script, and documents reruns. Run the focused pytest and expect failure because the guide does not yet exist.

- [ ] **Step 2: Write the example README**

Document:

1. Edit `config.py` once.
2. Run scripts 1, 2, either 3A or 3B, then 4.
3. Send `current_prompt.txt` to an LLM, review the YAML, and save it at `CELLTYPE_YAML`.
4. Run stage 5.
5. Rerun changed stages and all invalidated downstream stages.

Include POSIX and PowerShell command examples using the same filenames.

- [ ] **Step 3: Link the examples from the public docs**

Add a main README link to `examples/staged_workflow/README.md`. Add links from `docs/index.html` and `docs/QXYCell_overview.html` using the GitHub URL:

```text
https://github.com/jsonmad/QXYCell/tree/main/examples/staged_workflow
```

- [ ] **Step 4: Run focused tests and documentation checks**

Run the private structural pytest and expect all tests to pass. Parse all `docs/*.html`, validate repository-relative links, and confirm the public API names remain covered by the function reference and generated examples.

- [ ] **Step 5: Verify privacy and package boundaries**

Scan the new examples and changed docs for personal absolute paths, private dataset names, and email addresses. Build the wheel and source distribution, then assert neither archive contains `examples/`, `tests/`, or `docs/superpowers/`.

- [ ] **Step 6: Commit the guide and links**

Stage only the public example README and public documentation files; verify the local test remains ignored; commit with:

```bash
git commit -m "docs: link staged workflow examples"
```

---

### Task 3: Final integrated verification

**Files:**
- Verify only; no planned modifications.

**Interfaces:**
- Consumes: Tasks 1 and 2.
- Produces: evidence that the branch is documentation-integrated and build-safe.

- [ ] **Step 1: Run the private structural tests, Ruff, compilation, and import checks**

Expected: all commands exit zero and no numbered script runs its stage during import.

- [ ] **Step 2: Run HTML parsing, link validation, API coverage, and privacy scans**

Expected: all HTML parses, no missing links, no missing API names, and no private strings.

- [ ] **Step 3: Build and inspect distribution archives**

Expected: build succeeds and neither archive includes examples, tests, or internal design/plan files.

- [ ] **Step 4: Review the complete diff and branch status**

Expected: only intended public examples, documentation, and committed design/plan files are tracked; the private test remains ignored.
