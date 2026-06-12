import anndata as ad
import numpy as np
import pandas as pd

import qxycell as qxy


def test_plot_cell_boundaries_writes_png(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s1", "s1"],
            "celltype": ["A", "B", "A", "B"],
            "cell_polygon_wkt": [
                "POLYGON ((0 0, 8 0, 8 8, 0 8, 0 0))",
                "POLYGON ((10 0, 18 0, 18 8, 10 8, 10 0))",
                "POLYGON ((0 10, 8 10, 8 18, 0 18, 0 10))",
                "POLYGON ((10 10, 18 10, 18 18, 10 18, 10 10))",
            ],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)
    adata.obsm["spatial"] = np.array(
        [[4, 4], [14, 4], [4, 14], [14, 14]],
        dtype=float,
    )

    paths = qxy.plot_cell_boundaries(
        adata,
        sample_col="Sample",
        output_dir=tmp_path,
        label_celltypes="A",
        label_color="black",
        scale_bar=False,
        save_pdf=False,
        show=False,
        verbose=False,
    )

    assert len(paths["figures"]) == 1
    assert paths["figures"][0].name == "cell_boundaries_s1.png"
    assert paths["figures"][0].exists()
    assert paths["flip_y"] is True


def test_plot_cell_boundaries_accepts_raw_y_orientation(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1"],
            "celltype": ["A", "B"],
            "cell_polygon_wkt": [
                "POLYGON ((0 0, 8 0, 8 8, 0 8, 0 0))",
                "POLYGON ((10 0, 18 0, 18 8, 10 8, 10 0))",
            ],
        },
        index=["cell_0", "cell_1"],
    )
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)
    adata.obsm["spatial"] = np.array([[4, 4], [14, 4]], dtype=float)

    paths = qxy.plot_cell_boundaries(
        adata,
        sample_col="Sample",
        output_dir=tmp_path,
        flip_y=False,
        scale_bar=False,
        save_pdf=False,
        show=False,
        verbose=False,
    )

    assert paths["flip_y"] is False


def test_plot_cell_boundaries_defaults_to_polygon_bbox_centering(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1"],
            "celltype": ["A", "B"],
            "cell_polygon_wkt": [
                "POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))",
                "POLYGON ((100 0, 110 0, 110 10, 100 10, 100 0))",
            ],
        },
        index=["cell_0", "cell_1"],
    )
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)
    adata.obsm["spatial"] = np.array([[5, 5], [5, 5]], dtype=float)

    paths = qxy.plot_cell_boundaries(
        adata,
        sample_col="Sample",
        output_dir=tmp_path,
        scale_bar=False,
        save_pdf=False,
        show=False,
        verbose=False,
    )

    bounds = paths["sample_bounds"]["s1"]
    assert paths["center_method"] == "bbox"
    assert bounds["x_min"] == 0
    assert bounds["x_max"] == 110
    assert bounds["x_center"] == 55


def test_plot_cell_boundaries_underlay_false_centers_on_plotted_polygons(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1"],
            "celltype": ["A", "B"],
            "cell_polygon_wkt": [
                "POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))",
                "POLYGON ((1000 0, 1010 0, 1010 10, 1000 10, 1000 0))",
            ],
        },
        index=["cell_0", "cell_1"],
    )
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)
    adata.obsm["spatial"] = np.array([[5, 5], [1005, 5]], dtype=float)

    paths = qxy.plot_cell_boundaries(
        adata,
        sample_col="Sample",
        celltypes=["A"],
        underlay=False,
        output_dir=tmp_path,
        scale_bar=False,
        save_pdf=False,
        show=False,
        verbose=False,
    )

    bounds = paths["sample_bounds"]["s1"]
    assert bounds["x_min"] == 0
    assert bounds["x_max"] == 10
    assert bounds["x_center"] == 5
