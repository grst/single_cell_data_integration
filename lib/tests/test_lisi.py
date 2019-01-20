import numpy as np
import pandas as pd
from lib.lisi import inverse_simpson_index, lisi
import scanpy.api as sc
import pytest


def test_inversed_simpson_index():
    p_heter = np.array([1, 0, 0, 0, 0])

    for i in range(1, 20):
        assert inverse_simpson_index(np.ones(i)/i) == pytest.approx(i)
    assert inverse_simpson_index(p_heter) == 1


def test_lisi():
    X = np.random.rand(600, 10000)
    # init block matrix
    X[:300, :5000] += 0
    X[300:, 5000:] += 10
    no_mixture = ["A"] * 300 + ["B"] * 300
    perfect_mixture = ["A", "B"] * 300

    tdata = sc.AnnData(X=X)
    tdata.obs["no_mixture"] = no_mixture
    tdata.obs["perfect_mixture"] = perfect_mixture
    tdata.obs["no_mixture"] = tdata.obs["no_mixture"].astype("category")
    tdata.obs["perfect_mixture"] = tdata.obs["perfect_mixture"].astype("category")

    sc.pp.pca(tdata)

    lisi_no_mixture = lisi(tdata, "no_mixture")
    lisi_perfect_mixture = lisi(tdata, "perfect_mixture")

    assert np.median(lisi_no_mixture) == pytest.approx(1, 0.01)
    assert np.median(lisi_perfect_mixture) == pytest.approx(2, 0.01)
