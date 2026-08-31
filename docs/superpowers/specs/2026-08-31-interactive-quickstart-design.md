# Interactive Quickstart and Notebook Example Design

## Goal

Remove the redundant staged-workflow Python scripts and make the supported
quickstart paths clearly interactive: a persistent Python session or the
retained Jupyter notebook.

## Repository layout

- Delete every `.py` file from `staged_workflow/`, including `config.py`.
- Move `staged_workflow/QXYCell_staged_workflow.ipynb` to
  `examples/QXYCell_staged_workflow.ipynb` without changing its workflow.
- Move and rewrite `staged_workflow/README.md` as `examples/README.md`, focused
  only on opening and running the notebook.
- Remove the empty `staged_workflow/` directory.

## Quickstart documentation

Update the root README Quick start to explain that its code should be run in
one persistent interactive session. Tell readers to activate their QXYCell
environment, start `ipython` (preferred) or `python`, paste and execute each
stage in order, keep the session open so `adata` remains available, and pause
at the documented threshold and YAML review checkpoints. Link the notebook as
the cell-based alternative.

Update all checked-in Markdown and HTML references so none advertise or link
to the deleted numbered scripts. Existing documentation about the staged API,
checkpoint behavior, and rerun rules remains in scope and should not be
removed.

## Constraints

- Preserve the user's existing uncommitted edits in `README.md`.
- Do not change package behavior or public APIs.
- Do not modify notebook content except where a stale path or script reference
  requires correction.
- Keep changes direct and avoid adding new tooling or dependencies.

## Verification

- Confirm no `.py` files or stale `staged_workflow/` references remain.
- Confirm all local Markdown links affected by the move resolve.
- Parse the retained notebook as JSON and confirm its notebook structure is
  valid.
- Review the final Git diff to distinguish the user's pre-existing README edits
  from this task's changes.
