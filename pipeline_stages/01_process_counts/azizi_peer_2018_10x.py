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
import scanpy as sc
from anndata import AnnData
import os.path
from gene_identifiers import map_to_ensembl
import sys
sys.path.append("lib")
from scio import read_10x_mtx, concatenate, check_obs, check_var

DATASET = "azizi_peer_2018_10x"
OBS_DATA = "tables/datasets/{}_obs.tsv".format(DATASET)
MTX_BASENAME = "data/{}/{sample}_{patient}_TUMOR{replicate}"
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

obs = pd.read_csv(OBS_DATA, sep="\t")

# merge 10x
adatas = []
for i, row in obs.iterrows():
    filename = MTX_BASENAME.format(DATASET, sample=row['sample'], patient=row['patient'],
                                 replicate=row['replicate'])
    adata = read_10x_mtx(filename, var_names="gene_symbols")
    duplicated = adata.var_names.duplicated()
    print("Removing {} gene symbols because they are duplicated".format(sum(duplicated)))
    adata = adata[:, ~duplicated].copy()
    adata.obs["samples"] = row["sample"]
    adatas.append(adata)

adata = concatenate(adatas, merge_var_cols=["gene_ids"])

adata.obs = adata.obs.join(obs.set_index('sample'), on="samples", how="left")

adata.obs["dataset"] = DATASET

check_obs(adata)
check_var(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression="lzf")