import numpy as np
from numba import jit
import scanpy.api as sc


@jit(nopython=True)
def inverse_simpson_index(probs):
    """Compute the Inverse Simpson Index of an array or probabilities. """
    return 1/np.sum(np.array([x**2 for x in probs]))


@jit(nopython=True)
def _lisi(i, connectivities, annot_vec, levels):
    """Compute the locally inversed Simpson Index"""
    weights = connectivities[i, :]/np.sum(connectivities[i, :])
    probabs = []
    for batch in levels:
        mask = annot_vec == batch
        probabs.append(np.sum(weights[mask]))
    return inverse_simpson_index(probabs)


def lisi(adata, annot_col):
    """
    Compute the locally inversed Simpson Index (Harmony paper)

    Parameters:
        adata : AnnData
        annot_col : str
            index of the *categorical* column in `adata.obs`. E.g. `cell_type`
            or `batch`.
    """
    if 'lisi_connectivities' not in adata.uns:
        print("Computing connectivities using `sc.pp.neighbors`.")
        tmp_adata = sc.pp.neighbors(adata, n_neighbors=30, knn=False, method='gauss', copy=True)
        adata.uns['lisi_connectivities'] = tmp_adata.uns['neighbors']['connectivities']
        del tmp_adata
    annot_vec = adata.obs[annot_col].cat.codes.values
    levels = np.unique(annot_vec)
    return np.array([_lisi(i, adata.uns['lisi_connectivities'], annot_vec, levels)
                     for i in range(adata.shape[0])])

