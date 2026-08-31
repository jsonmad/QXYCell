# QXYCell interactive workflow example

The [QXYCell staged workflow notebook](QXYCell_staged_workflow.ipynb) runs the
core workflow one checkpoint at a time while keeping the active `adata` object
available between cells.

## Run the notebook

1. Activate the Python environment where QXYCell and Jupyter are installed.
2. From the repository root, run `jupyter lab` or `jupyter notebook`.
3. Open `examples/QXYCell_staged_workflow.ipynb`.
4. Edit the configuration cell, then run the cells in order.
5. In stage 3, choose exactly one threshold source: classifier JSON or a
   reviewed threshold table.
6. Pause after stage 4 for expert review of the generated cell-type YAML before
   applying it in stage 5.

Optional stage 2b removes cells in tissue or staining-artifact regions.
Optional stage 6 plots the assigned cell types.

## Rerun changed stages

Rerun the affected notebook cell and every dependent cell after it:

| Changed input | Rerun |
|---|---|
| Non-ignore annotation or cell GeoJSON | Stage 2, optionally stage 2b, the chosen stage 3 path, then stages 4 and 5 |
| Ignore annotation polygons | Stages 1, 2, and 2b, the chosen stage 3 path, then stages 4 and 5 |
| Classifier JSON thresholds | Stage 3 with classifier thresholds, then stages 4 and 5 |
| Reviewed threshold table | Stage 3 with the reviewed table, then stages 4 and 5 |
| Prompt context | Stage 4 |
| Cell-type YAML | Stage 5, then optionally stage 6 |
| Spatial plot settings | Stage 6 |

Stage 1 rebuilds the base checkpoint from measurements. Ignore-region changes
require this rebuild because stage 2b removes rows from the active checkpoint;
rerun stages 1 and 2 before removing cells again.
