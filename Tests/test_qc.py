import numpy as np
import pandas as pd
from anndata import AnnData


def test_low_gene_filter():

    X = np.random.poisson(1, (4, 100))

    obs = pd.DataFrame({
        "n_genes_by_counts": [50, 150, 300, 500],
        "pct_counts_mt": [5, 5, 5, 5]
    })

    adata = AnnData(X, obs=obs)

    filtered = adata[
        adata.obs["n_genes_by_counts"] >= 200
    ]

    assert filtered.n_obs == 2


def test_mitochondrial_filter():

    X = np.random.poisson(1, (4, 100))

    obs = pd.DataFrame({
        "n_genes_by_counts": [300, 300, 300, 300],
        "pct_counts_mt": [5, 10, 35, 40]
    })

    adata = AnnData(X, obs=obs)

    filtered = adata[
        adata.obs["pct_counts_mt"] < 20
    ]

    assert filtered.n_obs == 2


def test_scrublet_columns_exist():

    X = np.random.poisson(1, (3, 100))

    adata = AnnData(X)

    adata.obs["doublet_score"] = [0.1, 0.8, 0.2]
    adata.obs["predicted_doublet"] = [False, True, False]

    assert "doublet_score" in adata.obs.columns
    assert "predicted_doublet" in adata.obs.columns
    assert adata.obs["predicted_doublet"].sum() == 1