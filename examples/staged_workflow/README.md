# QXYCell staged workflow scripts

These scripts run the core QXYCell workflow one checkpoint at a time. Each
script has one job, so annotation GeoJSON, thresholds, prompt context, and
cell-type YAML can be revised without reimporting the measurements.

## Configure once

Activate the QXYCell environment, then edit [`config.py`](config.py):

- `PROJECT_DIR`: QuPath project folder containing the exported inputs.
- `OUTPUT_DIR`: fixed QXYCell output folder shared by every step.
- `THRESHOLD_TABLE`: reviewed threshold CSV/TSV used only by stage 3B.
- `CELLTYPE_YAML`: expert-reviewed cell-type logic applied by stage 5.
- `PIXEL_SIZE_UM`: verified square-pixel size used to scale GeoJSON coordinates.
- `CELLTYPE_CONTEXT`: tissue and project context added to the LLM prompt.

Use an absolute path for `PROJECT_DIR` and `OUTPUT_DIR` on your computer. The
fixed output folder lets every new Python process reload the same H5AD.

## First run

Run the scripts from this folder in numerical order.

### macOS or Linux

```bash
python 01_import_measurements.py
python 02_add_annotations.py

# Choose exactly one threshold source:
python 03a_threshold_from_classifiers.py
# OR: python 03b_threshold_from_table.py

python 04_generate_celltype_prompt.py
```

### Windows PowerShell

```powershell
python .\01_import_measurements.py
python .\02_add_annotations.py

# Choose exactly one threshold source:
python .\03a_threshold_from_classifiers.py
# OR: python .\03b_threshold_from_table.py

python .\04_generate_celltype_prompt.py
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

After cell typing, create one spatial PNG per image:

```bash
python 06_plot_spatial_celltypes.py
```

In PowerShell:

```powershell
python .\06_plot_spatial_celltypes.py
```

The script lists every `qxy.plot_spatial()` option explicitly so you can edit
the plot settings in one place. It keeps the QXYCell defaults except for
`show=False`, which prevents interactive plot windows during command-line or
remote runs. With the supplied options, PNG files are written to:

```text
OUTPUT_DIR/plots/spatial_celltypes/
```

## Rerun changed stages

Every successful script updates the same H5AD and exported tables. Rerun the
changed stage and any dependent stages:

| Changed input | Rerun |
|---|---|
| Annotation or cell GeoJSON | `02_add_annotations.py`, the chosen stage 3 script, then stages 4 and 5 |
| Classifier JSON thresholds | `03a_threshold_from_classifiers.py`, then stages 4 and 5 |
| Reviewed threshold table | `03b_threshold_from_table.py`, then stages 4 and 5 |
| Prompt context | `04_generate_celltype_prompt.py` |
| Cell-type YAML | `05_apply_celltypes.py`, then optionally `06_plot_spatial_celltypes.py` |
| Spatial plot settings | `06_plot_spatial_celltypes.py` |

Stage 1 rebuilds the base checkpoint from measurements. It is normally run
once unless the measurement files themselves change.

The two stage 3 scripts are deliberately independent. Classifier-only mode
ignores threshold tables; table-only mode uses the configured table and never
falls back to classifier JSON.
