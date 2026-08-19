# QXYCell staged workflow scripts

These scripts run the core QXYCell workflow one checkpoint at a time. Each
script has one job, so annotation GeoJSON, thresholds, prompt context, and
cell-type YAML can be revised without reimporting the measurements.

## Scripts or notebook

The [QXYCell staged workflow notebook](QXYCell_staged_workflow.ipynb) and the
numbered Python scripts implement the same staged workflow. Use the scripts
from a terminal or the notebook for an interactive, cell-by-cell experience.

In the notebook, edit the configuration cell, run stages 1 and 2, optionally
run stage 2b to remove cells in tissue or staining artifact regions, select
`"classifiers"` or `"table"` for stage 3, and run stage 4. Pause to review the
LLM-generated YAML with an expert, then run stage 5 to assign cell types and
optional stage 6 to plot them. The rerun rules below apply equally to notebook
cells and numbered scripts.

## Configure once

Activate the QXYCell environment, then edit [`config.py`](config.py):

- `PROJECT_DIR`: QuPath project folder containing the exported inputs.
- `OUTPUT_DIR`: fixed QXYCell output folder shared by every step.
- `THRESHOLD_TABLE`: reviewed threshold CSV/TSV used only by stage 3B.
- `CELLTYPE_YAML`: expert-reviewed cell-type logic applied by stage 5.
- `PIXEL_SIZE_UM`: verified square-pixel size used to scale GeoJSON coordinates.
- `IGNORE_ANNOTATION_TEXT`: case-insensitive text identifying exclusion annotations.
- `CELLTYPE_CONTEXT`: tissue and project context added to the LLM prompt.

Use an absolute path for `PROJECT_DIR` and `OUTPUT_DIR` on your computer. The
fixed output folder lets every new Python process reload the same H5AD.

## First run

Run the scripts from this folder in numerical order.

### macOS or Linux

```bash
python 01_import_cells.py
python 02_add_annotations.py

# Optional: remove cells in tissue/staining artifact annotations containing "ignore":
python 02b_remove_cells_by_annotations.py

# Choose exactly one threshold source:
python 03a_threshold_from_classifiers.py
# OR: python 03b_threshold_from_table.py

python 04_generate_celltype_prompt.py

# Pause here: review the generated YAML with an expert, then save it at CELLTYPE_YAML.
python 05_apply_celltypes.py
python 06_plot_spatial_celltypes.py  # Optional
```

### Windows PowerShell

```powershell
python .\01_import_cells.py
python .\02_add_annotations.py

# Optional: remove cells in tissue/staining artifact annotations containing "ignore":
python .\02b_remove_cells_by_annotations.py

# Choose exactly one threshold source:
python .\03a_threshold_from_classifiers.py
# OR: python .\03b_threshold_from_table.py

python .\04_generate_celltype_prompt.py

# Pause here: review the generated YAML with an expert, then save it at CELLTYPE_YAML.
python .\05_apply_celltypes.py
python .\06_plot_spatial_celltypes.py  # Optional
```

Stage 4 writes `OUTPUT_DIR/celltype/current_prompt.txt`. Copy that prompt into
an LLM, review and correct the returned YAML with an expert, and save the YAML
at the `CELLTYPE_YAML` path configured in `config.py`. Then apply it:

```bash
python 05_apply_celltypes.py
```

In PowerShell, the equivalent command is:

```powershell
python .\05_apply_celltypes.py
```

### Optional step 6: plot spatial cell types

After cell typing, create one spatial PNG per Sample label when available,
falling back to one PNG per QuPath image:

```bash
python 06_plot_spatial_celltypes.py
```

In PowerShell:

```powershell
python .\06_plot_spatial_celltypes.py
```

The script lists every `qxy.plot_spatial()` option explicitly so you can edit
the plot settings in one place. Its `sample_col=None` setting automatically
prefers usable `Sample` labels and otherwise uses `Image`; set it explicitly
to `"Sample"` or `"Image"` to force either grouping. It keeps the other
QXYCell defaults except for `show=False`, which prevents interactive plot
windows during command-line or remote runs. With the supplied options, PNG
files are written to:

```text
OUTPUT_DIR/plots/spatial_celltypes/
```

## Rerun changed stages

Every successful script updates the same H5AD and exported tables. Rerun the
changed stage and any dependent stages:

| Changed input | Rerun |
|---|---|
| Non-ignore annotation or cell GeoJSON | `02_add_annotations.py`, optionally stage 2b, the chosen stage 3 script, then stages 4 and 5 |
| Ignore annotation polygons | Stages 1, 2, and 2b, the chosen stage 3 script, then stages 4 and 5 |
| Classifier JSON thresholds | `03a_threshold_from_classifiers.py`, then stages 4 and 5 |
| Reviewed threshold table | `03b_threshold_from_table.py`, then stages 4 and 5 |
| Prompt context | `04_generate_celltype_prompt.py` |
| Cell-type YAML | `05_apply_celltypes.py`, then optionally `06_plot_spatial_celltypes.py` |
| Spatial plot settings | `06_plot_spatial_celltypes.py` |

Stage 1 rebuilds the base checkpoint from measurements. It is normally run
once unless the measurement files themselves change.

Stage 2b removes rows from the shared checkpoint. Running it repeatedly with
unchanged annotations is safe, but removed cells cannot be restored from that
checkpoint. If an ignore polygon is added, removed, or changed, rerun stages 1,
2, and 2b in order before continuing downstream.

Use stage 2b for annotations drawn around tissue folds, damaged tissue, debris,
edge artifacts, staining artifacts, or other regions that should not contribute
cells to the analysis.

The two stage 3 scripts are deliberately independent. Classifier-only mode
ignores threshold tables; table-only mode uses the configured table and never
falls back to classifier JSON. Classifier-only mode also writes the exact values
it applied to `OUTPUT_DIR/thresholds/classifier_thresholds.tsv`; rerunning stage
3A replaces that stable table.
