import numpy as np
from numba import jit
import scanpy.api as sc
import pandas as pd
import scipy
from multiprocessing import Pool


@jit(nopython=True)
def inverse_simpson_index(probs):
    """Compute the Inverse Simpson Index of an array or probabilities. """
    return 1/np.sum(np.array([x**2 for x in probs]))


def _lisi(i, connectivities, annot_vec, levels):
    """Compute the locally inversed Simpson Index"""
    weights = connectivities[i, :]/np.sum(connectivities[i, :])
    probabs = []
    for batch in levels:
        mask = annot_vec == batch
        w = weights[0, mask] if scipy.sparse.issparse(weights) else weights[mask]
        probabs.append(np.sum(w))
    return inverse_simpson_index(probabs)


def lisi_connectivities(adata, n_neighbors=30, type="gaussian"):
    """

    Parameters
    ----------
    adata : AnnData
    n_neighbors : int
        number of neighbors. Either as gaussian kernel width (`type = 'gaussian'`)
        or as number of fixed neighbors (`type = 'fixed'`). In the Harmony paper,
        they use a gaussian kernel with width of 30.
    type : str
        `gaussian` or `fixed`.

    Returns
    -------
    np.array (2d)
        adjacency matrix

    """
    if type == 'gaussian':
        tmp_adata = sc.pp.neighbors(adata, n_neighbors=n_neighbors, knn=False, method='gauss', copy=True)
        return tmp_adata.uns['neighbors']['connectivities']

    else:
        tmp_adata = sc.pp.neighbors(adata, n_neighbors=n_neighbors, copy=True)
        return tmp_adata.uns['neighbors']['connectivities']


def lisi(connectivities, labels, n_jobs=None):
    """
    Compute the locally inversed Simpson Index (Harmony paper)

    Parameters:
        connectivities : sparse CSR matrix
            adjacency matrix computed with `lisi_connectivities` or `sc.pp.neighbors`.
        labels : np.array
            numpy array containing a label for each cell.
    """
    p = Pool(n_jobs)
    annot_vec = pd.Series(labels, dtype="category").cat.codes.values
    levels = np.unique(annot_vec)
    return np.array(
        p.starmap(_lisi, [(i, connectivities, annot_vec, levels)
                          for i in range(len(labels))], chunksize=500)
    )
