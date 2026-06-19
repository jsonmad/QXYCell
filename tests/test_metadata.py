import anndata as ad
import numpy as np
import pandas as pd

import qxycell as qxy


def test_add_metadata_writes_h5ad_safe_partial_matches(tmp_path):
    adata = ad.AnnData(
        X=np.zeros((3, 1)),
        obs=pd.DataFrame(
            {"Sample": ["s1", "s2", "missing"]},
            index=["cell_0", "cell_1", "cell_2"],
        ),
    )
    metadata = pd.DataFrame(
        {
            "Sample": ["s1", "s2"],
            "group": ["Tumor", "NT"],
            "flag": [True, False],
            "age": [60, 70],
        },
    )

    summary = qxy.add_metadata(
        adata,
        metadata,
        sample_col="Sample",
        metadata_sample_col="Sample",
        output_dir=tmp_path,
        verbose=False,
    )

    assert summary["n_matched_samples"] == 2
    assert str(adata.obs["group"].dtype) == "category"
    assert str(adata.obs["flag"].dtype) == "boolean"
    assert adata.obs["age"].iloc[:2].tolist() == [60.0, 70.0]
    assert np.isnan(adata.obs["age"].iloc[2])

    adata.write_h5ad(tmp_path / "metadata_safe.h5ad")


def test_add_metadata_handles_missing_sample_keys(tmp_path):
    adata = ad.AnnData(
        X=np.zeros((3, 1)),
        obs=pd.DataFrame(
            {"Patient ID": pd.Series(["p1", pd.NA, "p2"], dtype="string").values},
            index=["cell_0", "cell_1", "cell_2"],
        ),
    )
    metadata = pd.DataFrame(
        {
            "Patient ID": ["p1", "p2"],
            "flag": [True, False],
        },
    )

    summary = qxy.add_metadata(
        adata,
        metadata,
        sample_col="Patient ID",
        metadata_sample_col="Patient ID",
        output_dir=tmp_path,
        verbose=False,
    )

    assert summary["missing_in_metadata"] == ["nan"]
    assert summary["n_matched_samples"] == 2
    assert str(adata.obs["flag"].dtype) == "boolean"
    assert pd.isna(adata.obs["flag"].iloc[1])

    adata.write_h5ad(tmp_path / "missing_key_metadata_safe.h5ad")
