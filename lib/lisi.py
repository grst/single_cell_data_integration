import numpy as np
import scanpy.api as sc
import pandas as pd
import scipy
from sklearn.preprocessing import normalize


def inverse_simpson_index(probs):
    """Compute the Inverse Simpson Index of an array or probabilities. """
    return 1/np.sum(np.power(probs, 2), axis=1)


def _lisi(connectivities, annot_vec, levels):
    # divide each row by sum of the row.
    is_sparse = scipy.sparse.issparse(connectivities)
    weights = normalize(connectivities, norm='l1', axis=1)

    probabs = []
    for batch in levels:
        mask = annot_vec == batch
        if is_sparse:
            probabs.append(np.array(np.sum(weights[:, mask], axis=1)))
        else:
            probabs.append(np.sum(weights[:, mask], axis=1)[:, np.newaxis])
    probabs = np.hstack(probabs)
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


def lisi(connectivities, labels):
    """
    Compute the locally inversed Simpson Index (Harmony paper)

    Parameters:
        connectivities : sparse CSR matrix
            adjacency matrix computed with `lisi_connectivities` or `sc.pp.neighbors`.
        labels : np.array
            numpy array containing a label for each cell.
    """
    annot_vec = pd.Series(labels, dtype="category").cat.codes.values
    levels = np.unique(annot_vec)
    return _lisi(connectivities, annot_vec, levels)
