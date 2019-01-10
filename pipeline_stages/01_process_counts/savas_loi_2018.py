# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.6.7
# ---

import pandas as pd
import scanpy.api as sc
from anndata import AnnData
import os.path
from gene_identifiers import map_to_ensembl
import sys
sys.path.append("lib")
from scio import read_10x_mtx, concatenate, check_obs, check_var

DATASET = "savas_loi_2018"
BASENAME = "data/{}/unmerged/{}"
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

samples = ["GSM3011853_tils20", "GSM3011854_tils32"]

# merge 10x
adatas = []
for i, sample in enumerate(samples):
    filename = BASENAME.format(DATASET, sample)
    adata = read_10x_mtx(filename, var_names="gene_symbols")
    adata.obs["samples"] = sample
    adata.obs["patient"] = str(i)
    adata.obs["origin"] = "tumor_primary"
    adata.obs["replicate"] = 1
    adata.obs["platform"] = "10x_3p_v2"
    adata.obs["tumor_type"] = "BRCA"
    duplicated = adata.var_names.duplicated()
    print("Removing {} gene symbols because they are duplicated".format(sum(duplicated)))
    adata = adata[:, ~duplicated].copy()
    adatas.append(adata)

adata = concatenate(adatas, merge_var_cols=["gene_ids"])

adata.obs["dataset"] = DATASET

check_obs(adata)
check_var(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
