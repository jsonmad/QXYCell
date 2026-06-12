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
