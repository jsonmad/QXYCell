# Basic step-by-step QXYCell script design

Date: 2026-08-09

## Purpose

Add a basic, directly editable Python script that performs a complete QXYCell analysis as explicit function calls instead of using `qxy.workflow()`. It should serve both as a runnable starting point and as concise documentation of the pipeline.

## Location

The script will be added as:

`scripts/basic_qxycell_step_by_step.py`

Because it is a single example script, no additional package module or abstraction will be introduced.

## User experience

The user edits a small configuration section at the top of the file and runs the script. Optional inputs use `None` or boolean switches, so the default path remains easy to understand.

Every analysis stage will begin with a numbered heading comment. Immediately beneath it, comments will explain:

- what the QXYCell function does;
- what it adds to or changes in `adata`;
- every parameter exposed by that function's current public signature, including its default and purpose;
- which values the example uses.

Optional stages will remain visible and runnable through explicit switches rather than being hidden in helper functions.

## Pipeline stages

1. Configure project paths and analysis options.
2. Call `qxy.check()` to validate the QuPath export and report inputs without creating AnnData.
3. Call `qxy.run()` to import measurements, spatial coordinates, annotations, polygons, and run metadata without implicitly thresholding or cell typing.
4. Call `qxy.threshold()` to apply the selected threshold table and create marker positivity columns.
5. Optionally call `qxy.add_metadata()` to attach sample-level metadata.
6. Call `qxy.remove_ignore()` to remove cells assigned to Ignore annotations.
7. Optionally call `qxy.celltype()` to apply ordered cell-type rules.
8. Call `qxy.qc()` to generate QC summaries and reports.
9. Call the standard plotting functions individually, including stacked bars, spatial plots, annotation polygons, marker positivity and intensity heatmaps, and cell boundaries where their required data are available.
10. Optionally call `qxy.cn_knn()`, `qxy.cn_kmeans()`, `qxy.cn_name()`, and CN plotting functions.
11. Call `qxy.save()` to write the final H5AD.

The ordering deliberately applies thresholds before cell typing, removes Ignore cells before QC and downstream analysis, and saves after all optional derived columns have been created.

## Error handling and optional inputs

- Required project paths will be represented by obvious placeholders that fail with a clear message if not edited.
- A failed `qxy.check()` will stop the example before importing data.
- Thresholding remains an explicit stage. The threshold-table path can be supplied or left as `None` to use QXYCell's normal threshold discovery behavior.
- Metadata and cell-type logic stages are skipped when their paths are `None`.
- Cell-type-dependent plots and CN analysis are skipped with a clear console message if `celltype` is unavailable.
- Annotation- or polygon-dependent plotting will rely on the plotting functions' validation and will be controlled by explicit switches.
- Output paths will reuse the active run directory stored in `adata.uns["qxycell"]` unless the example explicitly demonstrates an override.

## Scope

This change adds only the example script. It will not change QXYCell behavior, public function signatures, dependencies, or the existing notebook. It will not call the `qxy.workflow()` convenience wrapper.

## Verification

- Compare every documented option against the installed repository function signatures.
- Compile the script with Python to catch syntax errors.
- Run a lightweight static execution check that imports the script without launching a real analysis, if the script structure permits it.
- Run focused tests covering public API exports if necessary; no production behavior is being changed.
