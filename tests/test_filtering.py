import anndata as ad
import numpy as np
import pandas as pd
import pytest

import qxycell as qxy


def _adata_with_annotations():
    obs = pd.DataFrame(
        {
            "annotation__Ignore": [True, False, False, False],
            "annotation__Artifacts": [False, True, False, False],
            "annotation__Region": [False, False, True, False],
        },
        index=["cell_0", "cell_1", "cell_2", "cell_3"],
    )
    return ad.AnnData(X=np.zeros((4, 1)), obs=obs)


def _adata_with_labelled_annotations():
    obs = pd.DataFrame(
        {
            "annotation__IM": [True, False, False, False, False],
            "annotation__NT": [False, True, False, False, False],
            "annotation__Tonsil": [False, False, True, False, False],
            "annotation__Tumor": [False, False, False, True, False],
            "annotation__Ignore": [False, False, False, False, True],
        },
        index=["cell_0", "cell_1", "cell_2", "cell_3", "cell_4"],
    )
    adata = ad.AnnData(X=np.zeros((5, 1)), obs=obs)
    adata.uns["qxycell_annotation_labels"] = {
        "annotation__IM": "IM",
        "annotation__NT": "NT",
        "annotation__Tonsil": "Tonsil",
        "annotation__Tumor": "Tumor",
        "annotation__Ignore": "Ignore*",
    }
    return adata


def test_assign_annotations_assigns_labels_and_drops_columns():
    adata = _adata_with_labelled_annotations()

    summary = qxy.assign_annotations(
        adata,
        ["IM (22)", "NT (24)", "Tonsil (5)", "Tumor (26)"],
        target_col="Region",
        verbose=False,
    )

    assert adata.obs["Region"].astype(object).where(
        adata.obs["Region"].notna(), None
    ).tolist() == ["IM", "NT", "Tonsil", "Tumor", None]
    assert "annotation__IM" not in adata.obs.columns
    assert "annotation__Tumor" not in adata.obs.columns
    assert "annotation__Ignore" in adata.obs.columns
    assert "annotation__IM" not in adata.uns["qxycell_annotation_labels"]
    assert adata.uns["qxycell_annotation_labels"]["annotation__Ignore"] == "Ignore*"
    assert summary["target_col"] == "Region"
    assert summary["n_assigned_cells"] == 4
    assert summary["n_unassigned_cells"] == 1
    assert summary["dropped_annotation_columns"] is True


def test_assign_annotations_can_keep_source_columns():
    adata = _adata_with_labelled_annotations()

    qxy.assign_annotations(
        adata,
        "IM, Tumor",
        target_col="Region",
        drop=False,
        verbose=False,
    )

    assert adata.obs["Region"].astype(object).where(
        adata.obs["Region"].notna(), None
    ).tolist() == ["IM", None, None, "Tumor", None]
    assert "annotation__IM" in adata.obs.columns
    assert "annotation__Tumor" in adata.obs.columns


def test_assign_annotations_marks_overlapping_annotations_ambiguous():
    obs = pd.DataFrame(
        {
            "annotation__IM": [True, False],
            "annotation__Tumor": [True, True],
        },
        index=["cell_0", "cell_1"],
    )
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)

    with pytest.warns(UserWarning, match="more than one selected annotation"):
        summary = qxy.assign_annotations(
            adata,
            ["IM", "Tumor"],
            target_col="Region",
            verbose=False,
        )

    assert adata.obs["Region"].astype(str).tolist() == ["Ambiguous", "Tumor"]
    assert summary["n_conflicting_cells"] == 1


def test_assign_annotations_rejects_missing_labels():
    adata = _adata_with_labelled_annotations()

    with pytest.raises(KeyError, match="Missing expected annotation columns"):
        qxy.assign_annotations(
            adata,
            ["IM", "Missing"],
            target_col="Region",
            verbose=False,
        )


def test_assign_core_ids_from_measurements_prefers_tma_core_then_parent():
    obs = pd.DataFrame(
        {
            "TMA Core": ["G-1", "", pd.NA, "Unassigned"],
            "Parent": ["fallback_ignored", "H-2", "I-3", "J-4"],
        },
        index=["cell_0", "cell_1", "cell_2", "cell_3"],
    )
    adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)

    summary = qxy.assign_core_ids_from_measurements(adata, verbose=False)

    assert adata.obs["CoreID"].astype(object).where(
        adata.obs["CoreID"].notna(), None
    ).tolist() == ["G-1", "H-2", "I-3", "J-4"]
    assert summary["available_source_cols"] == ["TMA Core", "Parent"]
    assert summary["n_assigned_cells"] == 4
    assert summary["n_unassigned_cells"] == 0


def test_assign_core_ids_from_measurements_can_fill_only_missing_values():
    obs = pd.DataFrame(
        {
            "CoreID": ["existing", pd.NA, ""],
            "Parent": ["ignored", "H-2", "I-3"],
        },
        index=["cell_0", "cell_1", "cell_2"],
    )
    adata = ad.AnnData(X=np.zeros((3, 1)), obs=obs)

    qxy.assign_core_ids_from_measurements(
        adata,
        source_cols=("Parent",),
        fill_only_missing=True,
        verbose=False,
    )

    assert adata.obs["CoreID"].astype(object).where(
        adata.obs["CoreID"].notna(), None
    ).tolist() == ["existing", "H-2", "I-3"]


def test_assign_core_ids_from_measurements_rejects_missing_source_columns():
    adata = ad.AnnData(
        X=np.zeros((1, 1)),
        obs=pd.DataFrame({"Image": ["img"]}, index=["cell_0"]),
    )

    with pytest.raises(KeyError, match="None of the requested measurement core columns"):
        qxy.assign_core_ids_from_measurements(adata, verbose=False)


def test_remove_annotations_removes_cells_matching_text():
    adata = _adata_with_annotations()

    filtered = qxy.remove_annotations(
        adata,
        text="artifact",
        copy=True,
        verbose=False,
    )

    assert filtered.obs_names.tolist() == ["cell_0", "cell_2", "cell_3"]
    assert adata.obs_names.tolist() == ["cell_0", "cell_1", "cell_2", "cell_3"]


def test_remove_ignore_uses_ignore_text_wrapper():
    adata = _adata_with_annotations()

    qxy.remove_ignore(adata, verbose=False)

    assert adata.obs_names.tolist() == ["cell_1", "cell_2", "cell_3"]


def test_remove_annotations_rejects_empty_text():
    adata = _adata_with_annotations()

    with pytest.raises(ValueError, match="non-empty"):
        qxy.remove_annotations(adata, text="", verbose=False)
