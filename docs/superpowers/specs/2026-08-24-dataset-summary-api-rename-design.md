# Dataset summary API rename design

## Goal

Replace the misleading `qxy.qc()` API with `qxy.dataset_summary()`. The
function summarizes data already present in an AnnData object; it does not
validate image quality, segmentation accuracy, staining quality, thresholds,
batch effects, or spatial alignment.

## Public API

- Export `qxy.dataset_summary()` from `qxycell`.
- Remove `qxy.qc()` outright. Do not retain a compatibility alias.
- Move the implementation from `qxycell.qc` to `qxycell.summary`.
- Preserve the existing function arguments and returned dictionary structure,
  except for renamed output paths and metadata described below.

## Outputs

Write summary artifacts under `<run-folder>/dataset_summary/`:

- `dataset_summary.html`
- `dataset_overview.tsv`
- `cells_per_sample.tsv`, when the selected sample column exists
- `celltype_counts.tsv`, when the selected cell-type column exists
- `celltypes_per_sample.tsv`, when both columns exist
- `marker_positivity.tsv`, when marker-positive columns exist
- `marker_positivity_by_sample.tsv`, when marker-positive and sample columns exist
- `annotation_counts.tsv`, when annotation columns exist
- `ignore_cells.tsv`, when annotation columns containing `ignore` exist

Store output paths in `adata.uns["qxycell_dataset_summary"]`. Remove use of
`adata.uns["qxycell_qc"]`.

## Workflow integration

Rename the convenience-workflow option `make_qc` to
`make_dataset_summary`. When enabled, the workflow calls
`dataset_summary()`. Store the new option name in
`adata.uns["qxycell_workflow"]`; do not retain the old option.

## Documentation

- Rename the README heading and contents link from `QC` to `Dataset summary`.
- Describe the function as a descriptive summary, not a quality assessment.
- Update runnable examples and checked-in HTML documentation to use the new
  function, paths, metadata key, and terminology.

## Verification

- Add tests that exercise the new function and its renamed artifacts.
- Assert that `qxy.dataset_summary` is public and `qxy.qc` is absent.
- Update workflow tests for `make_dataset_summary`.
- Search tracked text files for obsolete public `qc` identifiers.
- Run the relevant test suite and documentation consistency checks.

## Scope

This change renames the existing summary capability without adding new quality
metrics or changing its calculations.
