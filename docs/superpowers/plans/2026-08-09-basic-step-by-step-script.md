# Basic Step-by-Step QXYCell Script Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a runnable Python example that performs all eleven QXYCell workflow stages through explicit public function calls and documents every option of every function it calls.

**Architecture:** Create one import-safe script with editable constants and a `main()` function. Keep each QXYCell operation visible in execution order, guard optional stages with booleans or nullable paths, and do not introduce package code or use `qxy.workflow()`.

**Tech Stack:** Python 3.11+, QXYCell public API, pathlib, pytest, Python AST/compile checks.

## Global Constraints

- Add only an example script and its focused test; do not change QXYCell public behavior, signatures, or dependencies.
- Cover all stages 1 through 11 from the approved design.
- Start every stage with a numbered `#` heading and explain what the function does and changes.
- List every parameter in the current public signature for every called QXYCell function, including defaults and meanings.
- Keep metadata, cell typing, plots, and CN analysis explicitly configurable.
- Do not call `qxy.workflow()`.

---

## File structure

- Create `scripts/basic_qxycell_step_by_step.py`: editable, executable pipeline example and inline option reference.
- Create `tests/test_basic_step_by_step_script.py`: verifies syntax, import safety, ordered API coverage, and absence of the workflow wrapper.

### Task 1: Add the executable step-by-step example

**Files:**
- Create: `tests/test_basic_step_by_step_script.py`
- Create: `scripts/basic_qxycell_step_by_step.py`

**Interfaces:**
- Consumes: public functions exported by `import qxycell as qxy`.
- Produces: `main() -> None`, with no analysis performed when the file is imported.

- [ ] **Step 1: Write the failing structure test**

Create `tests/test_basic_step_by_step_script.py` with a test that reads the example, compiles it, and verifies the required calls occur in pipeline order:

```python
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "basic_qxycell_step_by_step.py"


def test_basic_step_by_step_script_covers_full_pipeline():
    source = SCRIPT.read_text(encoding="utf-8")
    compile(source, str(SCRIPT), "exec")

    required_calls = [
        "qxy.check(",
        "qxy.run(",
        "qxy.threshold(",
        "qxy.add_metadata(",
        "qxy.remove_ignore(",
        "qxy.celltype(",
        "qxy.qc(",
        "qxy.plot_stacked_bar(",
        "qxy.plot_spatial(",
        "qxy.plot_annotation_polygons(",
        "qxy.plot_marker_positivity_heatmap(",
        "qxy.plot_marker_intensity_heatmap(",
        "qxy.plot_cell_boundaries(",
        "qxy.cn_knn(",
        "qxy.cn_kmeans(",
        "qxy.cn_name(",
        "qxy.plot_cn_heatmap(",
        "qxy.save(",
    ]
    positions = [source.index(call) for call in required_calls]
    assert positions == sorted(positions)
    assert "qxy.workflow(" not in source
    assert 'if __name__ == "__main__":' in source
```

- [ ] **Step 2: Run the test and verify it fails because the script is absent**

Run:

```bash
PYTHONPATH=src python3 -m pytest tests/test_basic_step_by_step_script.py -v
```

Expected: failure at `SCRIPT.read_text(...)` with `FileNotFoundError`.

- [ ] **Step 3: Create the script configuration and main guard**

Create `scripts/basic_qxycell_step_by_step.py` with:

```python
"""Run the complete QXYCell workflow as explicit, editable steps."""

from pathlib import Path

import qxycell as qxy


# 1. CONFIGURE PATHS AND ANALYSIS OPTIONS
PROJECT_DIR = Path("/replace/with/your/qupath/project")
OUTPUT_DIR: Path | None = None
THRESHOLD_FILE: Path | None = None
METADATA_FILE: Path | None = None
CELLTYPE_LOGIC: Path | None = None
SAMPLE_COL = "Sample"
METADATA_SAMPLE_COL: str | None = None
RUN_METADATA = False
RUN_CELLTYPING = False
RUN_PLOTS = True
RUN_CN_ANALYSIS = True
REMOVE_IGNORE_CELLS = True
K_NEIGHBORS = 12
N_CELLULAR_NEIGHBORHOODS = 12


def main() -> None:
    if not PROJECT_DIR.exists():
        raise FileNotFoundError(
            "Edit PROJECT_DIR at the top of this script before running it: "
            f"{PROJECT_DIR}"
        )

    # Stages 2-11 are inserted here as direct qxy calls.


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Add stages 2 through 8 with complete option comments**

Document and call these exact signatures, setting `qxy.run(..., apply_thresholds=False, celltype_logic=None)` so thresholding and cell typing remain visible separate steps:

```text
qxy.check(project_dir, output_dir=None, *, count_rows=False, threshold_file=None)
qxy.run(project_dir, output_dir=None, *, fail_on_check_error=True, pixel_size_um=0.28, threshold_file=None, apply_thresholds=False, celltype_logic=None, verbose=True)
qxy.threshold(adata, project_dir=None, *, threshold_file=None, output_dir=None, image_col="Image", verbose=True)
qxy.add_metadata(adata, metadata, *, sample_col="Image", metadata_sample_col=None, columns=None, prefix="", overwrite=False, output_dir=None, verbose=True)
qxy.remove_ignore(adata, *, annotation_prefix="annotation__", ignore_text="Ignore", copy=False, verbose=True)
qxy.celltype(adata, logic=None, *, celltype_column="celltype", unknown_label="Unknown", celltype_dir=None, verbose=True)
qxy.qc(adata, *, sample_col="Image", celltype_col="celltype", annotation_prefix="annotation__", marker_suffix="_pos", output_dir=None, verbose=True)
```

After `qxy.check()`, explicitly stop when `report.ok` is false. Because `remove_ignore(copy=False)` mutates AnnData in place, do not assign its return value. Skip metadata or cell typing with a printed explanation when the corresponding switch/path is absent.

- [ ] **Step 5: Add stage 9 plotting calls with complete option comments**

Call, in order, `plot_stacked_bar`, `plot_spatial`, `plot_annotation_polygons`, `plot_marker_positivity_heatmap`, `plot_marker_intensity_heatmap`, and `plot_cell_boundaries`. List every parameter from each current signature in the comments immediately before its call. Use the publication-oriented defaults already established by QXYCell: PDF for stacked bars and heatmaps; PNG for spatial, annotation-polygon QC, and cell-boundary raster plots. Pass `show=False` so batch execution does not open figures.

Only run cell-type plots when `celltype` exists. Only call cell-boundary plotting when `cell_polygon_wkt` exists. Keep annotation-polygons as explicit PNG-only QC output.

- [ ] **Step 6: Add stage 10 CN calls with complete option comments**

Document and call these signatures:

```text
qxy.cn_knn(adata, k=12, cell_type_col="celltype", image_col="Image")
qxy.cn_kmeans(adata, n_cn=12, key="cn", random_state=0, *, n_clusters=None)
qxy.cn_name(adata, key="cn", compaction=None, output_dir=None, verbose=True)
qxy.plot_spatial(adata, ..., category_col="cn", sample_col=SAMPLE_COL, show=False)
qxy.plot_cn_heatmap(adata, *, cn_col="cn", category_col=None, sample_col="Image", include_missing_samples=False, condition_col=None, normalize="sample", cluster_rows=True, cluster_cols=True, width="single", cmap=None, annotate=False, row_strip=False, dpi=600, save_png=False, save_pdf=True, save_svg=False, save_tiff=False, output_dir=None, show=True, verbose=True)
```

Skip the complete CN stage with a printed explanation if `celltype` is absent.

- [ ] **Step 7: Add stage 11 save call with complete option comments**

Document `qxy.save(adata, path=None, *, output_dir=None, verbose=True)`, call it once after every optional analysis stage, and print the returned H5AD path.

- [ ] **Step 8: Run the focused test**

Run:

```bash
PYTHONPATH=src python3 -m pytest tests/test_basic_step_by_step_script.py -v
```

Expected: `1 passed`.

- [ ] **Step 9: Verify every documented option against the source signatures**

Run an AST-based signature report over `src/qxycell/{checks,pipeline,metadata,filtering,celltyping,qc,plotting,neighbourhoods,io_utils}.py` and compare it line-by-line with the option blocks in the script. Expected: every parameter listed once in the relevant function's comment block, with no invented parameters.

- [ ] **Step 10: Run maintained tests**

Run:

```bash
PYTHONPATH=src python3 -m pytest tests -q
```

Expected: all maintained tests pass.

- [ ] **Step 11: Review and commit**

Run `git diff --check`, inspect `git diff -- scripts/basic_qxycell_step_by_step.py tests/test_basic_step_by_step_script.py`, then commit only the new script, its test, and this plan:

```bash
git add scripts/basic_qxycell_step_by_step.py tests/test_basic_step_by_step_script.py docs/superpowers/plans/2026-08-09-basic-step-by-step-script.md
git commit -m "Add step-by-step QXYCell workflow script"
```
