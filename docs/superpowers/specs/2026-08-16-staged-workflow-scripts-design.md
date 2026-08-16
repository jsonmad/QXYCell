# Public staged-workflow scripts design

## Goal

Provide small, public Python scripts that let a QXYCell user run each core stage
independently. The examples must demonstrate the rerunnable checkpoint workflow
without including private paths, test data, or research-specific settings.

## Location and files

Add the examples under `examples/staged_workflow/`:

- `README.md` explains setup, execution order, reruns, and the alternative
  threshold routes.
- `config.py` contains the user-editable paths and workflow settings.
- `01_import_measurements.py` creates the base AnnData checkpoint.
- `02_add_annotations.py` reloads that checkpoint and applies GeoJSON data.
- `03a_threshold_from_classifiers.py` applies classifier JSON thresholds only.
- `03b_threshold_from_table.py` applies a named threshold table only.
- `04_generate_celltype_prompt.py` creates the current LLM prompt.
- `05_apply_celltypes.py` applies the reviewed cell-type YAML.

The main README and documentation landing page will link to these examples.
The examples remain outside the Python package and PyPI source distribution.

## Configuration

`config.py` will use `pathlib.Path` and repository-safe placeholder paths. It
will define:

- `PROJECT_DIR`: the QuPath project folder.
- `OUTPUT_DIR`: one fixed QXYCell run folder shared by every script.
- `THRESHOLD_TABLE`: the reviewed CSV/TSV threshold table for option 3B.
- `CELLTYPE_YAML`: the reviewed cell-type YAML applied at stage 5.
- `PIXEL_SIZE_UM`: the verified square-pixel size used for GeoJSON scaling.
- `CELLTYPE_CONTEXT`: tissue and project context included in the LLM prompt.

Users edit this file once. A fixed `OUTPUT_DIR` makes separate Python processes
operate on the same H5AD without guessing which timestamped run is newest.

## Script behaviour

Each script will expose a `main()` function and invoke it only under
`if __name__ == "__main__"`. Stage 1 calls `qxy.import_measurements()` with the
configured project and output folders. Stages 2 through 5 call
`qxy.load(OUTPUT_DIR)` before invoking exactly one staged function.

The threshold scripts are explicit alternatives. Option 3A calls only
`qxy.threshold_from_classifiers()`. Option 3B calls only
`qxy.threshold_from_table()` with `THRESHOLD_TABLE`; neither script implements
fallback selection.

Stage 4 calls `qxy.celltype_prompt()` with `CELLTYPE_CONTEXT`. The user copies
the saved prompt to an LLM, reviews the returned YAML, and saves it at
`CELLTYPE_YAML`. Stage 5 passes that path explicitly to `qxy.celltype()`.

QXYCell's staged functions already checkpoint the H5AD and exported tables.
The example scripts will not add duplicate save calls or custom state handling.

## Rerun model

The README will tell users to run stages in numerical order for a first pass.
After measurement import, they may rerun a changed stage directly:

- New GeoJSON: rerun stage 2, then the selected threshold stage and later work.
- Changed thresholds: rerun either 3A or 3B, then stages 4 and 5.
- Changed prompt context: rerun stage 4.
- Changed cell-type YAML: rerun stage 5.

The examples rely on QXYCell's existing invalidation and replacement semantics.

## Error handling

The scripts will keep error handling direct. Missing paths, missing checkpoints,
conflicting classifiers, unresolved tables, or invalid YAML will surface through
QXYCell's existing descriptive exceptions. The examples will not catch and hide
these errors.

## Verification

Verification will include:

- A failing-first structural test for the expected example files and one-stage
  call mapping, followed by implementation and a passing test.
- Python compilation and Ruff checks for every example.
- Import-without-execution checks for each script.
- Documentation link validation.
- A scan for absolute personal paths, private dataset names, and email addresses.
- Existing documentation/API consistency checks.

No real QuPath dataset is bundled, so full data execution is outside this
example-only evidence boundary.
