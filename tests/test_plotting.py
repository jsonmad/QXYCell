import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import qxycell as qxy


def test_explicit_marker_heatmap_functions_are_public():
    assert callable(qxy.plot_marker_positivity_heatmap)
    assert callable(qxy.plot_marker_intensity_heatmap)


def test_plot_stacked_bar_can_save_pdf_only(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "celltype": ["A", "B", "A", "B"],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)

    paths = qxy.plot_stacked_bar(
        adata,
        sample_col="Sample",
        output_dir=tmp_path,
        save_png=False,
        save_pdf=True,
        show=False,
        verbose=False,
    )

    assert paths["png"] is None
    assert paths["pdf"].exists()
    assert paths["pdf"].suffix == ".pdf"


def test_plot_marker_heatmap_can_save_png_only(tmp_path):
    obs = pd.DataFrame(
        {
            "celltype": ["A", "A", "B", "B"],
            "CD3_pos": [True, False, True, True],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)
    adata.var_names = ["CD3"]

    paths = qxy.plot_marker_positivity_heatmap(
        adata,
        markers=["CD3"],
        cluster_rows=False,
        cluster_cols=False,
        save_png=True,
        save_pdf=False,
        save_svg=False,
        save_tiff=False,
        output_dir=tmp_path,
        show=False,
        verbose=False,
    )

    assert paths["positivity"][0].suffix == ".png"
    assert paths["positivity"][0].exists()
    assert paths["positivity"][1].suffix == ".csv"


def test_marker_intensity_heatmap_defaults_to_thresholded_markers(tmp_path):
    obs = pd.DataFrame(
        {
            "celltype": ["T", "T", "Other", "Other"],
            "CD3_pos": [True, True, False, False],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    adata = ad.AnnData(
        X=np.array([[10.0, 100.0], [8.0, 90.0], [1.0, 80.0], [2.0, 70.0]]),
        obs=obs,
    )
    adata.var_names = ["CD3", "UnthresholdedMarker"]

    paths = qxy.plot_marker_intensity_heatmap(
        adata,
        cluster_rows=False,
        cluster_cols=False,
        output_dir=tmp_path,
        show=False,
        verbose=False,
    )

    matrix = pd.read_csv(paths["intensity"][-1], index_col=0)
    assert matrix.columns.tolist() == ["CD3"]


def test_marker_intensity_heatmap_allows_explicit_unthresholded_markers(tmp_path):
    obs = pd.DataFrame(
        {"celltype": ["A", "B"]},
        index=["cell_0", "cell_1"],
    )
    adata = ad.AnnData(X=np.array([[1.0], [2.0]]), obs=obs)
    adata.var_names = ["RawMarker"]

    paths = qxy.plot_marker_intensity_heatmap(
        adata,
        markers=["RawMarker"],
        cluster_rows=False,
        cluster_cols=False,
        output_dir=tmp_path,
        show=False,
        verbose=False,
    )

    matrix = pd.read_csv(paths["intensity"][-1], index_col=0)
    assert matrix.columns.tolist() == ["RawMarker"]


def test_plot_annotation_polygons_falls_back_to_project_dir(tmp_path):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    (project_dir / "img.geojson").write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {
                            "objectType": "annotation",
                            "classification": {"name": "Region"},
                        },
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [[[0, 0], [20, 0], [20, 10], [0, 10], [0, 0]]],
                        },
                    },
                    {
                        "type": "Feature",
                        "properties": {"objectType": "cell"},
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [[[1, 1], [2, 1], [2, 2], [1, 2], [1, 1]]],
                        },
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    obs = pd.DataFrame(
        {"Image": ["img.tif", "img.tif", "img.tif"]},
        index=["cell_0", "cell_1", "cell_2"],
    )
    adata = ad.AnnData(X=np.zeros((3, 1)), obs=obs)
    adata.obsm["spatial"] = np.array([[2.0, 2.0], [5.0, 5.0], [8.0, 8.0]])
    adata.uns["qxycell"] = {
        "project_dir": str(project_dir),
        "output_dir": str(tmp_path / "run"),
        "pixel_size_um": 0.5,
    }

    result = qxy.plot_annotation_polygons(
        adata,
        output_dir=tmp_path / "plots",
        underlay_bins=32,
        show=False,
        verbose=False,
    )

    assert result["project_dir"] == project_dir.resolve()
    assert result["pixel_size_um"] == 0.5
    assert result["polygon_counts"] == {"img.tif": 1}
    assert result["cell_underlay"] is True
    assert result["fill"] is False
    assert result["underlay_bins"] == 32
    assert len(result["figures"]) == 1
    assert result["figures"][0].exists()


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


def test_plot_cell_boundaries_can_save_pdf_only(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1"],
            "celltype": ["A"],
            "cell_polygon_wkt": ["POLYGON ((0 0, 8 0, 8 8, 0 8, 0 0))"],
        },
        index=["cell_0"],
    )
    adata = ad.AnnData(X=np.zeros((1, 1)), obs=obs)
    adata.obsm["spatial"] = np.array([[4.0, 4.0]])

    paths = qxy.plot_cell_boundaries(
        adata,
        sample_col="Sample",
        output_dir=tmp_path,
        scale_bar=False,
        save_png=False,
        save_pdf=True,
        show=False,
        verbose=False,
    )

    assert len(paths["figures"]) == 1
    assert paths["figures"][0].suffix == ".pdf"
    assert paths["figures"][0].exists()


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


def test_plot_spatial_defaults_to_bbox_centering(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s1"],
            "celltype": ["A", "A", "B"],
        },
        index=["cell_0", "cell_1", "cell_2"],
    )
    adata = ad.AnnData(X=np.zeros((3, 1)), obs=obs)
    adata.obsm["spatial"] = np.array(
        [
            [0, 0],
            [10, 0],
            [100, 100],
        ],
        dtype=float,
    )

    paths = qxy.plot_spatial(
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
    assert bounds["x_center"] == 50
    assert bounds["y_center"] == 50


def test_plot_spatial_excludes_missing_sample_cells_by_default(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", np.nan],
            "celltype": ["A", "B", "A"],
        },
        index=["cell_0", "cell_1", "cell_missing"],
    )
    adata = ad.AnnData(X=np.zeros((3, 1)), obs=obs)
    adata.obsm["spatial"] = np.array([[0.0, 0.0], [10.0, 10.0], [1000.0, 1000.0]])

    paths = qxy.plot_spatial(
        adata,
        sample_col="Sample",
        output_dir=tmp_path,
        scale_bar=False,
        show=False,
        verbose=False,
    )

    assert list(paths["sample_bounds"]) == ["s1"]
    assert paths["sample_bounds"]["s1"]["x_max"] == 10.0
    assert paths["include_missing_samples"] is False
    assert paths["n_missing_sample_cells_excluded"] == 1
    assert len(paths["figures"]) == 1


def test_plot_spatial_accepts_different_underlay_adata(tmp_path):
    full_obs = pd.DataFrame(
        {"Sample": ["s1", "s1", "s1"], "celltype": ["A", "B", "Excluded"]},
        index=["cell_0", "cell_1", "cell_2"],
    )
    full = ad.AnnData(X=np.zeros((3, 1)), obs=full_obs)
    full.obsm["spatial"] = np.array([[0.0, 0.0], [10.0, 10.0], [100.0, 100.0]])
    plotted = full[:2].copy()
    plotted.obs["cn"] = ["N0", "N1"]

    paths = qxy.plot_spatial(
        plotted,
        underlay_adata=full,
        sample_col="Sample",
        category_col="cn",
        output_dir=tmp_path,
        scale_bar=False,
        show=False,
        verbose=False,
    )

    assert paths["underlay_adata_is_plot_adata"] is False
    assert paths["underlay_n_obs"] == 3
    assert paths["sample_bounds"]["s1"]["x_max"] == 100.0


def test_plot_spatial_can_save_pdf_only(tmp_path):
    obs = pd.DataFrame(
        {"Sample": ["s1", "s1"], "celltype": ["A", "B"]},
        index=["cell_0", "cell_1"],
    )
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)
    adata.obsm["spatial"] = np.array([[0.0, 0.0], [10.0, 10.0]])

    paths = qxy.plot_spatial(
        adata,
        sample_col="Sample",
        output_dir=tmp_path,
        scale_bar=False,
        save_png=False,
        save_pdf=True,
        show=False,
        verbose=False,
    )

    assert paths["save_png"] is False
    assert paths["save_pdf"] is True
    assert len(paths["figures"]) == 1
    assert paths["figures"][0].suffix == ".pdf"
    assert paths["figures"][0].exists()


def test_plot_cn_heatmap_accepts_category_col_alias(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "cn": ["A", "B", "A", "B"],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)

    paths = qxy.plot_cn_heatmap(
        adata,
        sample_col="Sample",
        category_col="cn",
        cluster_rows=False,
        cluster_cols=False,
        output_dir=tmp_path,
        show=False,
        verbose=False,
    )

    assert paths["sample"][0].name == "cn_heatmap_by_sample.pdf"
    assert paths["sample"][0].exists()
    assert paths["sample"][1].exists()
    assert paths["sample"][1].suffix == ".csv"


def test_plot_cn_heatmap_can_save_png_only(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "cn": ["A", "B", "A", "B"],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)

    paths = qxy.plot_cn_heatmap(
        adata,
        sample_col="Sample",
        cluster_rows=False,
        cluster_cols=False,
        save_png=True,
        save_pdf=False,
        save_svg=False,
        save_tiff=False,
        output_dir=tmp_path,
        show=False,
        verbose=False,
    )

    assert paths["sample"][0].suffix == ".png"
    assert paths["sample"][0].exists()
    assert paths["sample"][1].suffix == ".csv"


def test_plot_cn_heatmap_excludes_missing_samples_by_default(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", np.nan, "nan", "  ", "None", "<NA>"],
            "cn": ["A", "A", "B", "B", "A", "B"],
        },
        index=[f"cell_{i}" for i in range(6)],
    )
    adata = ad.AnnData(X=np.zeros((6, 1)), obs=obs)

    paths = qxy.plot_cn_heatmap(
        adata,
        sample_col="Sample",
        cluster_rows=False,
        cluster_cols=False,
        output_dir=tmp_path,
        show=False,
        verbose=False,
    )

    heatmap = pd.read_csv(paths["sample"][-1], index_col=0)
    assert heatmap.columns.tolist() == ["s1"]


def test_cn_heatmap_uses_vector_tiles(tmp_path, monkeypatch):
    import matplotlib.axes

    def reject_raster_heatmap(*args, **kwargs):
        raise AssertionError("Heatmap core must not be drawn with imshow().")

    monkeypatch.setattr(matplotlib.axes.Axes, "imshow", reject_raster_heatmap)
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "cn": ["A", "B", "A", "B"],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)

    paths = qxy.plot_cn_heatmap(
        adata,
        sample_col="Sample",
        cluster_rows=False,
        cluster_cols=False,
        output_dir=tmp_path,
        show=False,
        verbose=False,
    )

    assert paths["sample"][0].exists()


def test_plot_cn_heatmap_rejects_conflicting_category_alias(tmp_path):
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "cn": ["A", "B", "A", "B"],
            "other_cn": ["X", "Y", "X", "Y"],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)

    with pytest.raises(ValueError, match="Use either 'cn_col' or 'category_col'"):
        qxy.plot_cn_heatmap(
            adata,
            cn_col="other_cn",
            category_col="cn",
            output_dir=tmp_path,
            show=False,
            verbose=False,
        )
