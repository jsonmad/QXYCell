# QXYCell

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/qxycell-badge-dark.png">
    <img src="assets/qxycell-badge.png" alt="QXYCell badge" width="180">
  </picture>
</p>

QXYCell converts cell measurements and spatial assets from multiplex
immunofluorescence images processed in QuPath into an analysis-ready AnnData
`.h5ad` object for cell typing, plotting, and spatial analysis.

The resulting `.h5ad` object can be used with downstream tools such as
[Scanpy](https://scanpy.readthedocs.io/en/stable/),
[Squidpy](https://squidpy.readthedocs.io/en/stable/),
[scimap](https://scimap.xyz/), and
[SpatialData](https://spatialdata.scverse.org/en/stable/).

QXYCell is independent and is not affiliated with
[QuPath](https://qupath.github.io/)
([GPLv3](https://github.com/qupath/qupath/blob/main/LICENSE)).

## Workflow

![QXYCell workflow from multiplex tissue imaging through QuPath, QXYCell,
AnnData, spatial plots, and downstream analysis](docs/assets/qxycell_workflow.png)

Prepare the QuPath inputs first, then use the compact workflow below and open
the focused guides when more detail is needed.

## Installation

QXYCell requires Python 3.10 or newer.

```bash
git clone https://github.com/jsonmad/QXYCell.git
cd QXYCell
conda env create -f environment.yml
conda activate qxycell
```

Verify the installation:

```bash
python -c "import qxycell; print('QXYCell import OK')"
qxycell --help
```

### Updating

```bash
conda activate qxycell
git pull
conda env update -f environment.yml --prune
```

## Prepare data in QuPath

Before running QXYCell, follow the
[QuPath preparation guide](docs/qupath_preparation.md)
([PDF version](docs/QXYCell_QuPath_Preparation_Guide.pdf)). It covers:

- Sample, tissue-feature, and imaging-artifact annotations exported as GeoJSON
- Cell segmentation and cell measurements exported as `.tsv` or `.csv`
- Single-object classifiers for channel/marker thresholds (`.json`)
- Cell-boundary geometry exported as GeoJSON

The guide targets QuPath 0.7.0 and multiplex immunofluorescence data from any
acquisition platform. QXYCell supports square pixels only and uses a default
pixel size of 0.28 µm.

## Quick start

```python
import qxycell as qxy

project_dir = "/path/to/qupath_project"

# Optional preflight: inspect inputs without running the analysis.
report = qxy.check(project_dir)

# Stage 1: import cell measurements and create the AnnData checkpoint.
adata = qxy.import_cells(project_dir)

# Stage 2: add or refresh annotations and optional cell polygons.
# Annotation names containing "sample" define adata.obs["Sample"].
qxy.add_annotations(adata, pixel_size_um=0.28)

# Optional Stage 2b: choose the identifier string used in annotation names.
qxy.remove_cells(adata, remove_cells="ignore")
# qxy.remove_cells(adata, remove_cells="folded_tissue")

# Stage 3: choose one threshold source.
qxy.threshold_from_classifiers(adata)
# qxy.threshold_from_table(adata, "/path/to/thresholds.tsv")

# Stage 4: generate the prompt used to draft celltype_logic.yaml.
qxy.celltype_prompt(
    adata,
    context="Describe the tissue and expected populations",
)

# Pause for biology domain expert review, save the reviewed cell type YAML, then continue.
# Stage 5: assign cell types from cell type logic.
qxy.celltype(adata, "/path/to/celltype_logic.yaml")

# Stage 6: Sanity check - spatial plot assigned cell types.
qxy.plot_spatial(adata, category_col="celltype", show=True)
```

Classifier thresholding saves the applied values to
`thresholds/classifier_thresholds.tsv`. Table thresholding uses only the named
reviewed table.

Each successful stage updates the active `.h5ad` and refreshes
`tables/cells_obs.csv` and `tables/markers_var.csv`. If annotations are updated
after cells have been removed, rerun `adata = qxy.import_cells(project_dir)`
before refreshing annotations and removing cells again.

For runnable numbered scripts and an interactive notebook, see the
[staged workflow](staged_workflow/README.md).

## Documentation

| Guide | Use it for |
|---|---|
| [QuPath preparation](docs/qupath_preparation.md) ([PDF](docs/QXYCell_QuPath_Preparation_Guide.pdf)) | Preparing images, segmenting cells, measuring features, and exporting QuPath assets |
| [QuPath inputs, annotations, and thresholds](docs/qupath_inputs.md) | Input requirements, sample and removal annotations, pixel calibration, threshold sources, conflicts, and TMA cores |
| [Running the staged workflow](docs/running_qxycell.md) | Checkpoints, rerun rules, output folders, validation, and the optional single-call workflow |
| [Analysis workflows](docs/analysis.md) | Dataset summaries, metadata, cell typing, and cellular neighbourhoods |
| [Plotting](docs/plotting.md) | Spatial figures, cell boundaries, annotation polygons, bars, heatmaps, formats, and palettes |
| [AnnData structure and outputs](docs/anndata_and_outputs.md) | Stored fields, provenance, output files, and save/load behavior |
| [Staged scripts and notebook](staged_workflow/README.md) | Running one Python script per stage or working interactively |

Additional reference material:

- [QXYCell overview](docs/QXYCell_overview.html)
- [Function and command reference](docs/QXYCell_function_reference.html)
- [Synthetic function examples](docs/qxy_function_examples.html)
- [Documentation index](docs/index.html)

## Support and license

Report reproducible bugs and feature requests through
[GitHub Issues](https://github.com/jsonmad/QXYCell/issues). Do not attach
identifiable patient or research data, credentials, or other sensitive
material to a public issue.

Report suspected vulnerabilities privately as described in
[SECURITY.md](SECURITY.md).

QXYCell is released under the [MIT License](LICENSE).
