from __future__ import annotations

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

import qxycell as qxy


matplotlib.use("Agg")


def _thresholded_marker_adata() -> ad.AnnData:
    var_names = ["CD3", "CD3_2", "CD3_3", "CD3_4", "CD3_5", "CD3_6", "CD3_7", "CD3_8"]
    source_columns = [
        "Cell: CD3 - Cy5: Mean",
        "Cell: CD3 - Cy5: Median",
        "Cytoplasm: CD3 - Cy5: Mean",
        "Cytoplasm: CD3 - Cy5: Median",
        "Membrane: CD3 - Cy5: Mean",
        "Membrane: CD3 - Cy5: Median",
        "Nucleus: CD3 - Cy5: Mean",
        "Nucleus: CD3 - Cy5: Median",
    ]
    obs = pd.DataFrame(
        {
            "celltype": ["Low", "Low", "High", "High"],
            "CD3-nuc_pos": [0, 0, 1, 1],
        },
        index=[f"cell_{index}" for index in range(4)],
    )
    var = pd.DataFrame(
        {
            "source_measurement_column": source_columns,
            "threshold_marker_name": [""] * len(var_names),
            "positivity_column": [""] * len(var_names),
        },
        index=var_names,
    )
    var.loc["CD3_8", "threshold_marker_name"] = "CD3-nuc"
    var.loc["CD3_8", "positivity_column"] = "CD3-nuc_pos"
    X = np.zeros((4, len(var_names)))
    X[:, 0] = [10, 10, 1, 1]
    X[:, 7] = [1, 1, 10, 10]
    return ad.AnnData(X=X, obs=obs, var=var)


def _csv_path(result: dict[str, list], mode: str):
    return next(path for path in result[mode] if path.suffix == ".csv")


def test_default_heatmaps_use_thresholded_var_mapping(tmp_path):
    adata = _thresholded_marker_adata()

    positivity_result = qxy.plot_marker_positivity_heatmap(
        adata,
        markers=None,
        cluster_rows=False,
        cluster_cols=False,
        save_pdf=True,
        output_dir=tmp_path / "positivity",
        show=False,
        verbose=False,
    )
    intensity_result = qxy.plot_marker_intensity_heatmap(
        adata,
        markers=None,
        cluster_rows=False,
        cluster_cols=False,
        save_pdf=True,
        output_dir=tmp_path / "intensity",
        show=False,
        verbose=False,
    )

    positivity = pd.read_csv(_csv_path(positivity_result, "positivity"), index_col=0)
    intensity = pd.read_csv(_csv_path(intensity_result, "intensity"), index_col=0)

    assert intensity.columns.tolist() == ["CD3-nuc"]
    assert positivity.columns.tolist() == ["CD3-nuc"]
    assert intensity.loc["High", "CD3-nuc"] > intensity.loc["Low", "CD3-nuc"]
    assert positivity.loc["High", "CD3-nuc"] == 1.0
    assert positivity.loc["Low", "CD3-nuc"] == 0.0


def test_explicit_marker_keeps_existing_intensity_selection(tmp_path):
    adata = _thresholded_marker_adata()

    result = qxy.plot_marker_intensity_heatmap(
        adata,
        markers=["CD3"],
        cluster_rows=False,
        cluster_cols=False,
        save_pdf=True,
        output_dir=tmp_path,
        show=False,
        verbose=False,
    )

    intensity = pd.read_csv(_csv_path(result, "intensity"), index_col=0)
    assert intensity.columns.tolist() == ["CD3"]
    assert intensity.loc["High", "CD3"] < intensity.loc["Low", "CD3"]
