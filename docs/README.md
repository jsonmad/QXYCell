# QXYCell documentation

Prepare a QuPath multiplex-immunofluorescence project folder and convert its single-cell data into an analysis-ready AnnData object.

## Start with QuPath preparation

Verify image channels and square-pixel calibration, segment cells, and export the exact measurement and GeoJSON assets QXYCell expects. QXYCell defaults to 0.28 µm/pixel; supply the verified value when it differs.

[Read the preparation guide](qupath_preparation.md) · [Download the PDF](QXYCell_QuPath_Preparation_Guide.pdf)

## The core staged workflow

Build the analysis one checkpoint at a time. Every successful stage refreshes the active H5AD, `cells_obs.csv`, and `markers_var.csv`.

1.  **Measurements:** create the base AnnData.
2.  **Annotations:** add or replace GeoJSON-derived data.
3.  **Optional artifact regions:** remove cells inside `ignore` annotations drawn around tissue or staining artifacts.
4.  **Thresholds:** explicitly choose classifier JSON or one threshold table. Classifier mode saves the applied values as a reusable table.
5.  **Prompt:** generate the LLM prompt used to draft the YAML.
6.  **Cell types:** apply the reviewed YAML.
7.  **Optional plotting:** save spatial cell-type plots.

```python
import qxycell as qxy

adata = qxy.import_cells("/path/to/qupath_project")
qxy.add_annotations(adata, pixel_size_um=0.28)

# Optional: remove cells in tissue/staining artifact annotations
qxy.remove_cells(adata, remove_cells="ignore")

# Choose exactly one threshold source:
qxy.threshold_from_classifiers(adata)
# Saves/replaces thresholds/classifier_thresholds.tsv
# qxy.threshold_from_table(adata, "thresholds.tsv")

qxy.celltype_prompt(adata, context="Describe the tissue and expected populations")
# Pause here, review the generated YAML with an expert, then:
qxy.celltype(adata, "/path/to/celltype_logic.yaml")
qxy.plot_spatial(adata, category_col="celltype", show=False)
```

**Designed for iteration:** new GeoJSON, thresholds, prompts, or YAML replace that stage’s prior outputs and invalidate dependent downstream results. Because ignored-region removal deletes cells, changed ignore polygons require rerunning measurements, annotations, and ignored-region removal in order.

[Read the staged workflow guide](QXYCell_overview.md) · [See every stage parameter](QXYCell_function_reference.md)

## User guides

### Inputs and annotations

QuPath exports, sample and artifact annotations, thresholds, pixel calibration, and TMA cores.

[Prepare QXYCell inputs](qupath_inputs.md)

### Run QXYCell

Staged checkpoints, rerun rules, output folders, validation, and the workflow shortcut.

[Follow the staged workflow](running_qxycell.md)

### Sample metadata

Match experimental, clinical, and batch fields to images, samples, or TMA cores.

[Add sample metadata](metadata.md)

### Cell typing

Generate, review, apply, diagnose, and revise marker-positivity rules.

[Assign cell types](cell_typing.md)

### Cellular neighbourhoods

Build local composition profiles, cluster them, and review neighbourhood labels.

[Analyse neighbourhoods](cellular_neighbourhoods.md)

### Plotting

Spatial plots, cell boundaries, annotation polygons, bars, heatmaps, formats, and palettes.

[Create figures](plotting.md)

### AnnData and outputs

Stored fields, dataset summaries, provenance, output locations, flat-file exports, and save/load behavior.

[Inspect the data model](anndata_and_outputs.md)

## Reference documentation

### Overview

Follow the staged workflow, rerun rules, inputs, outputs, and AnnData model.

[Open overview](QXYCell_overview.md)

### Function reference

Python and command-line functions, parameters, and outputs.

[Open reference](QXYCell_function_reference.md)

### Function examples

Reproducible examples generated from a synthetic QuPath project folder.

[Open examples](qxy_function_examples.md)

### GitHub README

Project overview, installation, minimal quick start, and links to every guide.

[Open README](../README.md)
