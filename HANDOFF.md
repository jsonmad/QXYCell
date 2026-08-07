# QXYCell Project Handoff

Last updated: 2026-08-07  
Repository: `/Users/json/Desktop/iCloud/PROJECTS/260522_QXYCell`  
Branch: `codex/split-threshold-celltype-steps`  
Remote: `https://github.com/jsonmad/QXYCell.git`

## Purpose

QXYCell imports QuPath-exported cell measurement tables and GeoJSON annotations
into AnnData, applies marker thresholds and cell-type rules, performs cellular
neighbourhood analysis, and produces QC and publication-oriented plots.

The package is installed from `src/qxycell` and targets Python 3.10 or newer.
The primary local environment used for development is the Miniforge environment
named `qxycell`.

## Pre-cleanup Git state (historical)

This section records the state found before the approved August cleanup; it is
retained as provenance for the large generated-output removal.

The working tree is intentionally dirty and contains substantial work that has
not been committed. The current HEAD is:

```text
c812fa4 Document threshold output folder and stale cell type handling
```

Do not reset, clean, restore, or discard files without reviewing the working
tree. In particular:

- Many tracked files beneath old `260610_tests/outputs_*` directories are
  deleted. These deletions predate the latest work and should not be restored or
  committed without an explicit decision.
- Core source, tests, README, and generated documentation have local changes.
- `tests/test_measurements.py` and `tests/test_neighbourhoods.py` are new.
- `UWA_for_qxy_testing/`, `qupath_UWA/`, and
  `docs/_qxy_function_examples_build/` are untracked data/build directories.
- Several historical duplicate files have names containing ` (1)`. They are not
  normal import targets and should not be treated as the canonical source.

Before committing, inspect:

```bash
git status --short
git diff -- src/qxycell tests README.md docs
```

## Important behavior

### Import and saving

`qxy.run()` imports measurements, coordinates, annotations, and polygons into
AnnData. It automatically writes:

```text
<output_dir>/h5ad/qxycell.h5ad
```

There is currently no `save=False` option on `run()`.

Calling `qxy.save(adata)` later overwrites that default H5AD with the current
AnnData state. Pass a different `filename` to preserve both versions.

### Thresholding

Thresholding is separate from import:

```python
adata = qxy.run(project_dir, output_dir=output_dir)
qxy.threshold(
    adata,
    project_dir=project_dir,
    threshold_file=threshold_file,
    output_dir=output_dir,
)
```

`qxy.run()` defaults to `apply_thresholds=False`. Thresholding can be performed
inside `run()` or `workflow()` with:

```python
adata = qxy.run(
    project_dir,
    output_dir=output_dir,
    threshold_file=threshold_file,
    apply_thresholds=True,
)
```

`qxy.workflow()` now accepts both `threshold_file` and `apply_thresholds`.

When multiple QuPath classifiers define different thresholds for the same
measurement and image scope:

- the conflict is reported rather than resolved by filename order;
- `tables/classifier_conflicts.csv` records the candidates;
- generated threshold rows contain candidate provenance;
- per-image cells for conflicted rows are left blank;
- thresholding refuses to apply the row until every image threshold is numeric.

After reviewing the table, it can remain inside an output folder if supplied
explicitly:

```python
qxy.threshold(adata, threshold_file="/path/to/output/thresholds/thresholds.tsv")
```

Different image columns can contain different thresholds. Editing the table and
rerunning thresholding replaces marker positivity calls. Threshold-dependent
cell-type columns are archived with a `__stale_celltype` suffix; rerun
`qxy.celltype()` afterward.

Previous generated output folders are ignored during automatic measurement,
GeoJSON, and threshold-table discovery. The active output folder's thresholds
remain eligible. Classifier JSON discovery currently searches all JSON files
under the project directory, including generated output folders; this was
discussed and deliberately restored to the pre-change behavior at the user's
request.

### Core IDs and annotations

- Core IDs come only from the QuPath measurement-table column `TMA Core`.
- QXYCell does not derive CoreID from GeoJSON annotation labels.
- If `TMA Core` is absent, no CoreID is reported or created.
- Check and run reporting list every annotation name and explain which
  `adata.obs` column is created from it.
- Sample-labelled annotations populate `adata.obs["Sample"]`.
- Other annotation labels populate `annotation__<safe_label>` boolean columns.

### Metadata

Metadata may be applied in stages. Example:

```python
qxy.add_metadata(
    adata,
    r"Z:\Jason\260715_Debbie\metadata\CoreID.csv",
    sample_col="CoreID",
    metadata_sample_col="CoreID",
)

qxy.add_metadata(
    adata,
    r"Z:\Jason\260715_Debbie\metadata\PatientMetaData.csv",
    sample_col="Patient_ID",
    metadata_sample_col="Patient_ID",
)
```

The first table adds `Patient_ID`; the second uses it as its key.

### Ignore cells

Remove cells assigned to QuPath Ignore annotations with:

```python
qxy.remove_ignore(adata, verbose=True)
```

This changes the in-memory AnnData. Call `qxy.save()` to persist it.

### Plotting

- Plot functions support selecting PNG, PDF, SVG, and TIFF through explicit
  `save_*` arguments.
- Defaults produce one recommended format per plot type rather than always
  producing both PNG and PDF.
- Annotation polygon QC plots default to `fill=False` and do not produce PDFs
  by default.
- Annotation polygon plots can use a low-resolution cell-density underlay.
- `plot_spatial()` supports `underlay_adata` for plotting one AnnData selection
  over cell locations from another.
- `plot_spatial()` and `plot_cn_heatmap()` exclude missing sample labels by
  default.
- Marker positivity and marker intensity heatmaps have separate public
  functions.
- PDF/SVG heatmap data grids use vector `pcolormesh` tiles rather than raster
  `imshow`, keeping the heatmap sharp when enlarged.

### Cellular neighbourhoods

Typical sequence:

```python
qxy.cn_knn(adata, k=12)
qxy.cn_kmeans(adata, n_cn=12)
qxy.cn_name(adata)
```

CN labels are Windows-safe. Generated labels replace `/` and `\` with ` + `,
replace Windows-reserved characters, collapse `..`, strip trailing periods and
spaces, and protect reserved device names.

## Windows UTF-8 fix

QuPath measurement TSV/CSV files are read explicitly with `utf-8-sig`.
This supports ordinary UTF-8 and UTF-8 files with a BOM, while preventing
`Centroid X µm` and `Centroid Y µm` from becoming `Âµ` on Windows.

The encoding is applied consistently to:

- measurement header validation;
- measurement row scans for images and TMA Core values;
- full pandas measurement imports.

Invalid bytes retain the prior `errors="replace"` behavior in CSV-reader paths.

## Main modules

```text
src/qxycell/checks.py          project validation and reports
src/qxycell/classifiers.py     classifier and threshold parsing
src/qxycell/measurements.py    measurement discovery and validation
src/qxycell/pipeline.py        run and threshold application
src/qxycell/celltyping.py      cell-type rules
src/qxycell/filtering.py       Ignore, sample, and CoreID handling
src/qxycell/geojson.py         GeoJSON discovery and geometry import
src/qxycell/metadata.py        sample/patient metadata joins
src/qxycell/neighbourhoods.py  KNN, KMeans, and CN naming
src/qxycell/plotting.py        spatial, bar, polygon, and heatmap plots
src/qxycell/qc.py              QC tables and reports
src/qxycell/io_utils.py        H5AD load/save helpers
src/qxycell/workflow.py        convenience workflow wrapper
```

## Verification

The last maintained-suite run completed successfully:

```text
75 passed in 3.45s
```

Command:

```bash
PYTHONPYCACHEPREFIX=/tmp/qxy-pycache \
MPLCONFIGDIR=/tmp/qxy-mplconfig \
LOKY_MAX_CPU_COUNT=4 \
conda run -n qxycell pytest -q \
  tests/test_celltyping.py \
  tests/test_filtering.py \
  tests/test_io_utils.py \
  tests/test_measurements.py \
  tests/test_metadata.py \
  tests/test_neighbourhoods.py \
  tests/test_output_paths.py \
  tests/test_plotting.py \
  tests/test_threshold_tables.py \
  tests/test_tma_geojson.py
```

Do not use unconstrained `pytest -q` as the authoritative suite yet. Repository
discovery also collects historical/manual files such as `test_heatmaps.py` and
duplicate `tests/* (1).py` files, which rely on external data or removed APIs.

## Windows deployment

A portable source snapshot was previously created under `dist/` with a
Miniforge installation guide. Because additional changes have been made since
that snapshot, regenerate it before treating it as the current Windows release.

For a shared OVD installation, prefer a non-editable, administrator-maintained
environment:

```powershell
conda create --prefix "C:\CondaCommon\envs\qxycell" python=3.12 pip -y
conda activate "C:\CondaCommon\envs\qxycell"
python -m pip install .
```

Grant users read/execute permission and restrict environment modification to
maintainers. Versioned shared environments are safer for updates that must not
interrupt active users.

## Suggested next steps

1. Review the large dirty working tree and decide which old generated-output
   deletions belong in version control.
2. Remove or formally exclude historical duplicate ` (1)` source/test files.
3. Add pytest configuration so manual external-data scripts are not collected.
4. Review and commit the source, tests, README, and generated docs as coherent
   changes.
5. Regenerate and validate the Windows source bundle after the commit.
6. Consider adding `save=False` to `qxy.run()` only if explicitly requested.

