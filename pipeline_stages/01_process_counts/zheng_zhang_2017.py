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
from gene_identifiers import map_to_ensembl, collapse_gene_symbols
import sys
sys.path.append("lib")
from scio import check_obs, check_var

DATASET = "zheng_zhang_2017"
COUNT_FILE = "data/{}/GSE98638_HCC.TCell.S5063.count.txt.gz".format(DATASET)
OBS_FILE = "data/{}/patient_table.tsv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)


# first position of SampleType
origin_map = {
    "P": "blood_peripheral",
    "N": "normal_adjacent",
    "T": "tumor_primary",
    "J": "tumor_edge"
}

# second and third
cell_type_map = {
    "TC": "T cell",
    "TH": "T cell CD4+ non-regulatory",
    "TR": "T cell regulatory (Treg)",
    "TY": "T cell CD25 int."
}

obs = pd.read_csv(OBS_FILE, sep="\t")
obs.columns = [c.strip() for c in obs.columns]

obs = obs.assign(origin = obs["sampleType"].apply(lambda x: origin_map[x[0]]))\
         .assign(cell_type = obs["sampleType"].apply(lambda x: cell_type_map[x[1:]]))\
         .assign(tumor_type = "LIHC")\
         .assign(replicate = 1)\
         .assign(platform = "smartseq2")\
         .rename({"Patient": "patient"}, axis="columns")

obs = obs.assign(samples = obs[["patient", "origin", "replicate"]].apply(
    lambda x: "_".join([str(k) for k in x]), axis=1))

obs = obs.set_index("UniqueCell_ID")

raw_counts = pd.read_csv(COUNT_FILE, sep="\t")

var = raw_counts[["geneID", "symbol"]].set_index("symbol")

adata = AnnData(raw_counts.iloc[:, 2:].values.transpose(), obs, var)
adata = adata[:, ~pd.isna(adata.var_names)]
# adata = map_to_ensembl(adata)

adata.obs["dataset"] = DATASET

check_obs(adata)
check_var(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
