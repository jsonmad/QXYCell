"""CN (cell neighbourhood) analysis using scikit-learn KNN and KMeans."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def cn_knn(
    adata: "AnnData",
    k: int = 12,
    cell_type_col: str = "celltype",
    image_col: str = "Image",
) -> "AnnData":
    """Build a CN composition profile using k-nearest neighbours.

    For each cell, finds its *k* nearest spatial neighbours within the same image
    (to avoid neighbours bleeding across tissue sections), then computes the
    frequency of each cell type among those neighbours. The result is a per-cell
    composition vector stored in ``adata.obsm["cn_profile"]``.

    Requires scikit-learn (``pip install scikit-learn``).

    Parameters
    ----------
    adata:
        AnnData object. Must have ``adata.obsm["spatial"]``, ``adata.obs[image_col]``,
        and ``adata.obs[cell_type_col]``.
    k:
        Number of nearest neighbours per cell (default 12). Self is excluded.
    cell_type_col:
        Column in ``adata.obs`` containing cell type labels (default ``"celltype"``,
        as assigned by ``qxy.celltype()``).
    image_col:
        Column in ``adata.obs`` used to separate cells by image (default ``"Image"``).
        KNN graphs are built independently per image.

    Returns
    -------
    AnnData
        The same object with:

        - ``adata.obsm["cn_profile"]`` — float32 array of shape
          ``(n_cells, n_cell_types)``. Each row sums to 1 (or 0 for isolated cells).
        - ``adata.uns["cn"]`` — dict with run parameters and the ordered list of
          cell types.
    """
    try:
        from sklearn.neighbors import NearestNeighbors
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for cn_knn. "
            "Install it with: pip install scikit-learn"
        ) from exc

    import numpy as np

    # --- validate ---
    if cell_type_col not in adata.obs.columns:
        raise ValueError(
            f"cell_type_col '{cell_type_col}' not found in adata.obs. "
            f"Available columns: {adata.obs.columns.tolist()}"
        )
    if image_col not in adata.obs.columns:
        raise ValueError(
            f"image_col '{image_col}' not found in adata.obs. "
            f"Available columns: {adata.obs.columns.tolist()}"
        )
    if "spatial" not in adata.obsm:
        raise ValueError(
            "adata.obsm['spatial'] not found. Run qxy.run() first to populate spatial coordinates."
        )

    # --- setup ---
    cell_types = sorted(adata.obs[cell_type_col].dropna().unique().tolist())
    ct_index = {ct: i for i, ct in enumerate(cell_types)}
    n_cells = adata.n_obs
    n_types = len(cell_types)

    profile = np.zeros((n_cells, n_types), dtype=np.float32)

    obs_ct = adata.obs[cell_type_col].values
    obs_img = adata.obs[image_col].values

    spatial_raw = adata.obsm["spatial"]
    # obsm values may be ndarray or DataFrame
    spatial = spatial_raw.values if hasattr(spatial_raw, "values") else np.asarray(spatial_raw)

    global_indices = np.arange(n_cells)
    images = np.unique(obs_img)

    # --- per-image KNN ---
    for image in images:
        mask = obs_img == image
        img_global = global_indices[mask]
        img_spatial = spatial[mask]
        img_ct = obs_ct[mask]
        n_img = len(img_global)

        if n_img <= 1:
            continue  # cannot form any neighbours

        # k+1 because kneighbors returns self as the first hit
        n_neighbors = min(k + 1, n_img)

        nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean", algorithm="auto")
        nn.fit(img_spatial)
        _, indices = nn.kneighbors(img_spatial)

        # indices[:, 0] is always self — start from 1
        for local_i, neighbour_local in enumerate(indices[:, 1:]):
            global_i = img_global[local_i]
            for ni in neighbour_local:
                ct = img_ct[ni]
                if ct in ct_index:
                    profile[global_i, ct_index[ct]] += 1

    # --- normalise rows to frequencies ---
    row_sums = profile.sum(axis=1, keepdims=True)
    nonzero = row_sums[:, 0] > 0
    profile[nonzero] = profile[nonzero] / row_sums[nonzero]

    adata.obsm["cn_profile"] = profile
    adata.uns["cn"] = {
        "k": k,
        "cell_type_col": cell_type_col,
        "image_col": image_col,
        "cell_types": cell_types,
    }

    print(
        f"cn_knn: {n_cells:,} cells | k={k} | {n_types} cell types | "
        f"{len(images)} image(s) → adata.obsm['cn_profile']"
    )
    return adata


def cn_kmeans(
    adata: "AnnData",
    n_cn: int = 12,
    key: str = "cn",
    random_state: int = 0,
    *,
    n_clusters: int | None = None,  # legacy alias for n_cn
) -> "AnnData":
    """Cluster cells into CNs using MiniBatchKMeans.

    Clusters ``adata.obsm["cn_profile"]`` (built by :func:`cn_knn`) into
    *n_cn* CNs and stores the labels as a ``pd.Categorical`` in
    ``adata.obs[key]``.

    Requires scikit-learn (``pip install scikit-learn``).

    Parameters
    ----------
    adata:
        AnnData object with ``adata.obsm["cn_profile"]`` present.
    n_cn:
        Number of CNs (default 12). Also accepted as ``n_clusters`` for
        backwards compatibility.
    key:
        Column name written to ``adata.obs`` (default ``"cn"``).
    random_state:
        Random seed passed to MiniBatchKMeans for reproducibility (default 0).

    Returns
    -------
    AnnData
        The same object with:

        - ``adata.obs[key]`` — ``pd.Categorical`` of labels ``"N0"`` … ``"N{n-1}"``.
        - ``adata.uns["cn"]`` — updated with clustering parameters.
    """
    try:
        from sklearn.cluster import MiniBatchKMeans
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for cn_kmeans. "
            "Install it with: pip install scikit-learn"
        ) from exc

    import pandas as pd

    if n_clusters is not None:
        n_cn = n_clusters  # honour legacy kwarg

    if "cn_profile" not in adata.obsm:
        raise ValueError(
            "adata.obsm['cn_profile'] not found. Run qxy.cn_knn(adata) first."
        )

    profile = adata.obsm["cn_profile"]

    kmeans = MiniBatchKMeans(
        n_clusters=n_cn,
        random_state=random_state,
        n_init="auto",
    )
    labels = kmeans.fit_predict(profile)

    adata.obs[key] = pd.Categorical([f"N{label}" for label in labels])

    if "cn" not in adata.uns:
        adata.uns["cn"] = {}
    adata.uns["cn"].update(
        {
            "n_cn": n_cn,
            "key": key,
            "method": "MiniBatchKMeans",
            "random_state": random_state,
        }
    )

    print(f"cn_kmeans: {n_cn} CNs → adata.obs['{key}']")
    return adata


def cn_name(
    adata: "AnnData",
    key: str = "cn",
    compaction: dict[str, str] | None = None,
    output_dir=None,
    verbose: bool = True,
) -> "object":
    """Apply deterministic short labels to CN clusters from their composition profile.

    Labels are derived from ``adata.obsm["cn_profile"]`` using the following rules
    applied in order for each CN:

    +-------------------------------------+-----------------------+-------------------------+
    | Condition                           | Label form            | Example                 |
    +=====================================+=======================+=========================+
    | Top cell type ≥ 50 %                | ``<top> hi``          | ``Tumor hi``            |
    +-------------------------------------+-----------------------+-------------------------+
    | Top cell type ≥ 35 %                | ``<top> mix``         | ``Macrophage mix``      |
    +-------------------------------------+-----------------------+-------------------------+
    | Top two cell types together ≥ 55 %  | ``<top1>/<top2>``     | ``CD8 T/Macrophage``    |
    +-------------------------------------+-----------------------+-------------------------+
    | Otherwise                           | ``<top1>/<top2> mix`` | ``CD8 T/NK mix``        |
    +-------------------------------------+-----------------------+-------------------------+

    ``adata.obs[key]`` is updated in place: cluster IDs (``N0``, ``N1``, …) are
    replaced with descriptive labels. The original ID→label mapping is preserved in
    ``adata.uns["cn"]["label_map"]``.

    Parameters
    ----------
    adata:
        AnnData object with ``adata.obsm["cn_profile"]`` and ``adata.obs[key]``
        present (i.e. after running :func:`cn_knn` and :func:`cn_kmeans`).
    key:
        Column in ``adata.obs`` containing CN cluster IDs (default ``"cn"``).
    compaction:
        Optional dict mapping original cell type names to shorter display labels,
        e.g. ``{"CD8+PD1+LAG3+": "PD1 LAG3 CD8"}``. Applied before naming.
        Unmatched cell types are used as-is.
    output_dir:
        Folder to write ``cn_labels.csv``. Defaults to the QXYCell output folder
        stored in ``adata.uns``.
    verbose:
        Print a summary table (default True).

    Returns
    -------
    pandas.DataFrame
        Label table with columns ``CN_ID``, ``CN_Label``, ``Top_Contributors``.
        Also saved to ``<output_dir>/cn/cn_labels.csv``.
    """
    import numpy as np
    import pandas as pd

    from qxycell.paths import resolve_output_dir

    # --- validate ---
    if key not in adata.obs.columns:
        raise ValueError(
            f"'{key}' not found in adata.obs. Run qxy.cn_kmeans(adata) first."
        )
    if "cn_profile" not in adata.obsm:
        raise ValueError(
            "adata.obsm['cn_profile'] not found. Run qxy.cn_knn(adata) first."
        )
    if "cn" not in adata.uns or "cell_types" not in adata.uns["cn"]:
        raise ValueError(
            "adata.uns['cn']['cell_types'] not found. Run qxy.cn_knn(adata) first."
        )

    cell_types = adata.uns["cn"]["cell_types"]
    compaction = compaction or {}
    display_types = [compaction.get(ct, ct) for ct in cell_types]

    profile = adata.obsm["cn_profile"]
    if hasattr(profile, "values"):
        profile = profile.values
    profile = np.asarray(profile, dtype=np.float32)

    obs_cn = adata.obs[key].astype(str).values
    cn_ids = list(adata.obs[key].cat.categories) if hasattr(adata.obs[key], "cat") else sorted(set(obs_cn))

    # --- per-CN mean composition ---
    composition: dict[str, "np.ndarray"] = {}
    for cn_id in cn_ids:
        mask = obs_cn == str(cn_id)
        composition[cn_id] = profile[mask].mean(axis=0) if mask.any() else np.zeros(len(cell_types))

    # --- apply naming rules ---
    def _label(comp: "np.ndarray") -> str:
        sorted_idx = comp.argsort()[::-1]
        top1_i = sorted_idx[0]
        top2_i = sorted_idx[1] if len(sorted_idx) > 1 else None
        top1_f = float(comp[top1_i])
        top1_n = display_types[top1_i]

        if top1_f >= 0.50:
            return f"{top1_n} hi"
        if top1_f >= 0.35:
            return f"{top1_n} mix"
        if top2_i is not None and (top1_f + float(comp[top2_i])) >= 0.55:
            return f"{top1_n}/{display_types[top2_i]}"
        top2_n = display_types[top2_i] if top2_i is not None else ""
        return f"{top1_n}/{top2_n} mix" if top2_n else f"{top1_n} mix"

    raw_labels = {cn_id: _label(composition[cn_id]) for cn_id in cn_ids}

    # --- disambiguate duplicate labels ---
    seen: dict[str, int] = {}
    label_map: dict[str, str] = {}
    for cn_id in cn_ids:
        lbl = raw_labels[cn_id]
        if lbl in seen:
            seen[lbl] += 1
            label_map[cn_id] = f"{lbl} ({seen[lbl]})"
        else:
            seen[lbl] = 1
            label_map[cn_id] = lbl

    duplicates = [lbl for lbl, count in seen.items() if count > 1]
    if duplicates and verbose:
        print(f"cn_name: duplicate labels disambiguated: {duplicates}")

    # --- build output table ---
    rows = []
    for cn_id in cn_ids:
        comp = composition[cn_id]
        sorted_idx = comp.argsort()[::-1]
        top_contributors = "; ".join(
            f"{display_types[i]} {comp[i]:.0%}"
            for i in sorted_idx[:3]
            if comp[i] > 0.01
        )
        rows.append({
            "CN_ID": cn_id,
            "CN_Label": label_map[cn_id],
            "Top_Contributors": top_contributors,
        })
    label_table = pd.DataFrame(rows)

    # --- update adata ---
    new_categories = [label_map[cn_id] for cn_id in cn_ids]
    adata.obs[key] = pd.Categorical(
        adata.obs[key].astype(str).map(label_map),
        categories=new_categories,
    )
    adata.uns["cn"]["label_map"] = label_map
    adata.uns["cn"]["compaction"] = compaction

    # --- save table ---
    out_dir = (
        resolve_output_dir(adata=adata) / "cn"
        if output_dir is None
        else __import__("pathlib").Path(output_dir)
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "cn_labels.csv"
    label_table.to_csv(csv_path, index=False)

    if verbose:
        print(f"cn_name: {len(cn_ids)} CNs labelled → adata.obs['{key}']")
        print(label_table.to_string(index=False))
        print(f"  saved → {csv_path}")

    return label_table
