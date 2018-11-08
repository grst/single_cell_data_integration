import scanpy.api as sc
import pandas as pd
import numpy as np

def read_10x_mtx(basename, var_names="gene_symbols"):
    """
    Read 10x mtx files with basename.

    Usually, 10x mtx files are stored in
    a single folder per sample which can
    be natively read by scanpy.

    However, publicly available data often looks like this
      * GSM1234_genex.tsv
      * GSM1234_matrix.mtx
      * GSM1234_barcodes.tsv

    This function allows to read this 'basename' file format.


    Args:
        basename: the sample path before `_{matrix,barcodes/genes}`
        var_names: {'gene_symbols', 'gene_ids'}. Gene IDs = ENSG.

    """

    mtx_file, barcode_file, genes_file = ["{}_{}".format(basename, part)
            for part in ("matrix.mtx", "barcodes.tsv", "genes.tsv")]

    # read mtx
    adata = sc.read(mtx_file).T

    # read barcodes
    barcodes = pd.read_csv(barcode_file, sep="\t", header=None)
    adata.obs_names = barcodes[0]

    # read genes
    genes = pd.read_csv(genes_file, sep="\t", header=None)
    if var_names == "gene_symbols":
        adata.var_names = genes[1]
        adata.var['gene_ids'] = genes[0].values
    else:
        adata.var_names = genes[0]
        adata.var['gene_symbols'] = genes[1].values

    return adata


def concatenate(adatas, merge_var_cols=None, **kwargs):
    """
    Extension of scanpy's native `concatenate` funcion.
    Allows to merge columns in `var` of the same name with same
    contents into a single one instead of
    generating col-1, col-2, ...

    The columns to merge need to be specified explicitly.

    Args:
        adatas: list of `AnnData` objects
        merge_var_cols: columns of adata.var to merge into one. Default: None
        **kwargs: will be passed to `scanpy.concatenate`.

    """

    cols = {}
    merge_var_cols = [] if merge_var_cols is None else merge_var_cols

    # store a copy of the columns and check they are all identical
    for col_name in merge_var_cols:
        cols[col_name] = adatas[0].var[col_name]
        for adata in adatas[1:]:
            assert adata.shape[1] == adatas[0].shape[1], "shapes incompatible"
            assert np.all(adata.var[col_name] == cols[col_name]),\
                    "columns are not identical in all adata objects: {}".format(col_name)

    # remove columns before merging
    for col_name in merge_var_cols:
        for adata in adatas:
            adata.var.drop(col_name, axis="columns", inplace=True)

    adata = adatas[0].concatenate(adatas[1:], **kwargs)

    # re-add columns
    for col_name in merge_var_cols:
        adata.var[col_name] = cols[col_name]

    return adata


def check_obs(adata):
    """
    Check that the columns in obs follow
    the naming conventions.
    """
    obs = adata.obs
    mandatory_cols = ["sample", "patient", "origin",
            "replicate", "platform", "tumor_type"]
    for col in mandatory_cols:
        assert col in obs.columns, "{} is a mandatory column".format(col)
        assert np.sum(obs[col].isnull()) == 0, "NAs in column {}".format(col)

    # check sample columns
    sample_count = obs.groupby(["patient", "origin", "replicate"])['sample'].nunique()

    assert np.all(sample_count == 1), \
            "sample must be unique for each patient, origin and replicate"

    # check CV
    assert np.all(obs["origin"].isin(["tumor_primary", "normal_adjacent",
        "tumor_edge", "blood_peripheral"])), "invalid word in column origin"

    assert np.all(obs["platform"].isin(["10x_3p", "10x_3p_v2", "10x_5p",
        "indrop_v2", "smartseq2"])), "invalid word in column platform"


    tcga = ["LAML", "ACC", "BLCA", "LGG", "BRCA",
        "CESC", "CHOL", "LCML", "COAD", "CNTL",
        "ESCA", "FPPP", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC",
        "DLBC", "MESO", "MISC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM",
        "STAD", "TGCT", "THYM", "THCA", "UCS", "UCEC", "UVM"]
    other = ["PBMC", "NSCLC"]
    assert np.all(obs["tumor_type"].isin(tcga + other)), \
                "invalid word in column tumor_type"


def check_var(adata):
    """
    Check that the adata.var follows the conventions of thie project
    """
    assert adata.var_names.is_unique, "Var names must be unique"

    assert np.all(x[:4] == ENSG for x in adata.var_names), \
            "var names must be ensemble gene identifiers"

    assert "gene_symbols" in adata.var.columns, \
            "var must contain a gene_symbols column"
