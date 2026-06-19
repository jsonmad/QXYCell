import anndata as ad
import numpy as np
import pandas as pd
import pytest

import qxycell as qxy


def test_save_clears_obs_index_name_collision(tmp_path):
    obs = pd.DataFrame({"quxy_cell_id": ["cell_0", "cell_1"]})
    obs.index = obs["quxy_cell_id"].astype(str)
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)
    adata.obs["quxy_cell_id"] = adata.obs["quxy_cell_id"].astype("string").astype("category")

    path = qxy.save(adata, tmp_path / "collision_safe.h5ad", verbose=False)

    assert path.exists()
    assert adata.obs.index.name is None


def test_save_rejects_obs_index_name_collision_when_values_differ(tmp_path):
    obs = pd.DataFrame({"quxy_cell_id": ["different_0", "different_1"]})
    obs.index = pd.Index(["cell_0", "cell_1"], name="quxy_cell_id")
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)

    with pytest.raises(ValueError, match="index values do not match"):
        qxy.save(adata, tmp_path / "collision_unsafe.h5ad", verbose=False)
