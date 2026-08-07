import anndata as ad
import numpy as np
import pandas as pd
import pytest

import qxycell as qxy


def test_cn_name_generates_windows_safe_labels(tmp_path):
    adata = ad.AnnData(
        X=np.zeros((4, 1)),
        obs=pd.DataFrame(
            {"cn": pd.Categorical(["N0", "N0", "N1", "N1"])},
            index=[f"cell_{i}" for i in range(4)],
        ),
    )
    adata.obsm["cn_profile"] = np.array(
        [
            [0.8, 0.2],
            [0.7, 0.3],
            [0.2, 0.8],
            [0.3, 0.7],
        ],
        dtype=np.float32,
    )
    adata.uns["cn"] = {
        "cell_types": ["CD8/T..:*?", r"Macrophage\Tumour"],
    }

    labels = qxy.cn_name(adata, output_dir=tmp_path, verbose=False)

    forbidden = set('<>:"/\\|?*')
    for label in labels["CN_Label"]:
        assert not forbidden.intersection(label)
        assert ".." not in label
        assert not label.endswith((" ", "."))
    assert labels["CN_Label"].str.contains(" + ", regex=False).all()
    assert set(adata.obs["cn"].astype(str)) == set(labels["CN_Label"])
    assert (tmp_path / "cn_labels.csv").exists()


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("CD8+", "CD8 +"),
        ("Something.", "Something"),
        ("Something  ", "Something"),
        ("Something . ", "Something"),
    ],
)
def test_windows_safe_cn_label_has_no_trailing_dot_or_space(raw, expected):
    from qxycell.neighbourhoods import _windows_safe_cn_label

    assert _windows_safe_cn_label(raw) == expected
