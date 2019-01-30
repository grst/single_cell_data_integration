import numpy as np
import pandas as pd
from lib.lisi import inverse_simpson_index, lisi, lisi_connectivities
import scanpy.api as sc
import pytest
import numpy.testing as npt
import scipy.sparse


def test_inversed_simpson_index():
    p_heter = np.array([[1, 0, 0, 0, 0]])

    assert inverse_simpson_index(p_heter) == 1

    for i in range(1, 20):
        assert inverse_simpson_index(np.ones((1, i))/i)[0] == pytest.approx(i)

    # test with full array
    npt.assert_almost_equal(
        inverse_simpson_index(np.array([[1, 0, 0, 0], [0.5, 0.5, 0, 0]])),
        np.array([1, 2]),
        decimal=2
    )


def test_lisi_gaussian():
    X = np.random.rand(600, 10000)
    # init block matrix
    X[:300, :5000] += 0
    X[300:, 5000:] += 10
    no_mixture = ["A"] * 300 + ["B"] * 300
    perfect_mixture = ["A", "B"] * 300

    adata = sc.AnnData(X=X)
    sc.pp.pca(adata)

    connectivities = lisi_connectivities(adata)

    lisi_no_mixture = lisi(connectivities, no_mixture)
    lisi_perfect_mixture = lisi(connectivities, perfect_mixture)

    assert np.median(lisi_no_mixture) == pytest.approx(1, 0.01)
    assert np.median(lisi_perfect_mixture) == pytest.approx(2, 0.01)


def test_lisi_gaussian_sparse():
    X = np.random.rand(600, 10000)
    # init block matrix
    X[:300, :5000] += 0
    X[300:, 5000:] += 10
    no_mixture = ["A"] * 300 + ["B"] * 300
    perfect_mixture = ["A", "B"] * 300

    adata = sc.AnnData(X=scipy.sparse.csr_matrix(X))
    sc.pp.pca(adata)

    connectivities = lisi_connectivities(adata)

    lisi_no_mixture = lisi(connectivities, no_mixture)
    lisi_perfect_mixture = lisi(connectivities, perfect_mixture)

    assert np.median(lisi_no_mixture) == pytest.approx(1, 0.01)
    assert np.median(lisi_perfect_mixture) == pytest.approx(2, 0.01)


def test_lisi_fixed():
    X = np.random.rand(600, 10000)
    # init block matrix
    X[:300, :5000] += 0
    X[300:, 5000:] += 10
    no_mixture = ["A"] * 300 + ["B"] * 300
    perfect_mixture = ["A", "B"] * 300

    adata = sc.AnnData(X=X)
    sc.pp.pca(adata)

    connectivities = lisi_connectivities(adata, n_neighbors=100, type="fixed")

    lisi_no_mixture = lisi(connectivities, no_mixture)
    lisi_perfect_mixture = lisi(connectivities, perfect_mixture)

    assert np.median(lisi_no_mixture) == pytest.approx(1, 0.01)
    # with the fixed width (depending on the size of the fixed neighborhood)
    # slightly smaller values than two are normal.
    assert np.median(lisi_perfect_mixture) == pytest.approx(2, 0.05)


# def test_lisi_parallel():
#     X = np.random.rand(10000, 500)
#     # init block matrix
#     X[:5000, :250] += 0
#     X[5000:, 250:] += 10
#     perfect_mixture = ["A", "B"] * 15000
#
#     adata = sc.AnnData(X=X)
#     sc.pp.pca(adata)
#
#     connectivities = lisi_connectivities(adata)
#
#     print("run lisi with 8 jobs")
#     lisi_perfect_mixture = lisi(connectivities, perfect_mixture, n_jobs=8)

