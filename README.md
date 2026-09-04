# QXYCell

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/qxycell-badge-dark.png">
    <img src="docs/assets/qxycell-badge.png" alt="QXYCell badge" width="96">
  </picture>
</p>

A Python bridge from a QuPath project folder to AnnData. QXYCell converts cell measurements and spatial data from multiplex
immunofluorescence images processed in QuPath
([QuPath](https://qupath.github.io/);
[GPLv3](https://github.com/qupath/qupath/blob/main/LICENSE))
into an AnnData object for cell typing, visualization, and spatial analysis. The resulting `.h5ad` object can be used with downstream tools.

## Workflow

![QXYCell workflow from multiplex tissue imaging through QuPath, QXYCell,
AnnData, spatial plots, and downstream analysis](docs/assets/qxycell_workflow.png)

## Installation

QXYCell requires Python 3.10 or newer.

```bash
# Install in new environment
git clone https://github.com/jsonmad/QXYCell.git
cd QXYCell
conda env create -f environment.yml
conda activate qxycell
```
```bash
# Verify the installation
python -c "import qxycell; print('QXYCell import OK')"
qxycell --help
```
```bash
# update if required
conda activate qxycell
git pull
conda env update -f environment.yml --prune
```

## Prepare data in QuPath

Before running QXYCell, follow the
[QuPath preparation guide](docs/qupath_preparation.md). It covers:

- Sample, tissue-feature, and imaging-artifact annotations exported as GeoJSON
- Cell segmentation and cell measurements exported as `.tsv` or `.csv`
- Single-object classifiers for channel/marker thresholds (`.json`)
- Cell-boundary geometry exported as GeoJSON

## Quick start

To begin, define the path to the QuPath project directory and an output directory for the QXYCell results. The output directory can be anywhere but should not sit inside of the QuPath project directory.

Run this quickstart in an interactive Python session or convert to a jupyter notebook. Each data-processing stage checkpoints the current `adata` object to the active output folder. After activating the QXYCell environment, start `ipython` (recommended) or `python` from a terminal, then paste and execute each stage below in order. Pause where noted to review thresholds and cell-type YAML.

```console
ipython
# Or, if IPython is not installed:
python
```

At the interactive prompt, paste the Python code below, one stage at a time.


```python
import qxycell as qxy

project_dir = "/path/to/qupath_project"
output_dir = "/path/to/outputs/run_1"

# Optional preflight: inspect inputs without running the analysis.
report = qxy.check(project_dir)

# Stage 1: import cell measurements and create the AnnData checkpoint.
adata = qxy.import_cells(project_dir, output_dir=output_dir)

# Stage 2: add or refresh annotations and optional cell polygons.
# Annotation names containing "sample" define adata.obs["Sample"].
qxy.add_annotations(adata, pixel_size_um=0.28)

# Optional Stage 2b: choose the identifier string used in annotation names.
qxy.remove_cells(adata, remove_cells="ignore")
# qxy.remove_cells(adata, remove_cells="folded_tissue")

# Stage 3: choose either 3A or 3B. Do not run both.
# Classifier thresholding saves the applied values to thresholds/classifier_thresholds.tsv. Table thresholding uses only the named reviewed table.

# Stage 3A: apply thresholds from QuPath object-classifier JSON files.
qxy.threshold_from_classifiers(adata)

# Stage 3B: generate, review, and apply a threshold table instead.
threshold_table = qxy.generate_threshold_table(
    project_dir,
    output_dir=output_dir,
    )

# Pause here to review and fill every per-image threshold in the generated TSV.
# qxy.threshold_from_table(adata, threshold_table)

# Stage 4: generate the prompt used to draft celltype_logic.yaml.
qxy.celltype_prompt(
    adata,
    context="Describe the tissue and expected populations",
    )

# Pause for biology domain expert review, save the reviewed cell type YAML, then continue.
# Stage 5: assign cell types from cell type logic.
qxy.celltype(adata, "/path/to/celltype_logic.yaml")

# Stage 6: sanity check the assigned cell types spatially.
qxy.plot_spatial(adata, category_col="celltype", show=True)

# Stage 7: plot marker positivity and intensity by assigned cell type.
qxy.plot_marker_positivity_heatmap(
    adata,
    category_col="celltype",
    show=True,
    )

qxy.plot_marker_intensity_heatmap(
    adata,
    category_col="celltype",
    show=True,
    )

# Stage 8: compare marker positivity and intensity between samples.
qxy.plot_marker_positivity_heatmap(
    adata,
    category_col="Sample",
    show=True,
    )

qxy.plot_marker_intensity_heatmap(
    adata,
    category_col="Sample",
    show=True,
    )
```

You can exit `exit()` after any stage finishes successfully. To restart, activate the same environment, start a new interactive session, recreate the path variables, and load the `.h5ad` from the `output_dir` or the exact `.h5ad` path.

```python
# restarting a session
import qxycell as qxy

project_dir = "/path/to/qupath_project"
output_dir = "/path/to/outputs/run_1"
adata = qxy.load(output_dir)

# or exact path to .h5ad
adata = qxy.load("/path/to/outputs/run_1/h5ad/qxycell.h5ad")
```


Classifier thresholding saves the applied values to `thresholds/classifier_thresholds.tsv`. Table thresholding uses only the named
reviewed table.

Each successful stage updates the active `.h5ad` and refreshes
`tables/cells_obs.csv` and `tables/markers_var.csv`. If annotations are updated
after cells have been removed, rerun
`adata = qxy.import_cells(project_dir, output_dir=output_dir)` before refreshing
annotations and removing cells again.

## Documentation

| Guide | Use it for |
|---|---|
| [QuPath preparation](docs/qupath_preparation.md) | Preparing images, segmenting cells, measuring features, and exporting QuPath assets |
| [QuPath inputs, annotations, and thresholds](docs/qupath_inputs.md) | Input requirements, sample and removal annotations, pixel calibration, threshold sources, conflicts, and TMA cores |
| [Running the staged workflow](docs/running_qxycell.md) | Checkpoints, rerun rules, output folders, validation, and the optional single-call workflow |
| [Sample metadata](docs/metadata.md) | Matching experimental, clinical, and batch metadata to images, samples, or TMA cores |
| [Cell typing](docs/cell_typing.md) | Prompt generation, reviewed YAML rules, assignment diagnostics, validation, and reruns |
| [Cellular neighbourhoods](docs/cellular_neighbourhoods.md) | Local composition profiles, clustering, naming, parameter review, and neighbourhood plots |
| [Plotting](docs/plotting.md) | Spatial figures, cell boundaries, annotation polygons, bars, heatmaps, formats, and palettes |
| [AnnData structure and outputs](docs/anndata_and_outputs.md) | Stored fields, dataset summaries, provenance, output files, and save/load behavior |

Additional reference material:

- [QXYCell overview](docs/QXYCell_overview.md)
- [Function and command reference](docs/QXYCell_function_reference.md)
- [Synthetic function examples](docs/qxy_function_examples.md)
- [Documentation index](docs/README.md)

## Support and license

Report reproducible bugs and feature requests through
[GitHub Issues](https://github.com/jsonmad/QXYCell/issues).
Report suspected vulnerabilities privately as described in
[SECURITY.md](SECURITY.md).
QXYCell is released under the [MIT License](LICENSE).
