import pandas as pd
import scanpy.api as sc
from anndata import AnnData
import os.path
from gene_identifiers import map_to_ensembl
import sys
sys.path.append("lib")
from scio import read_10x_mtx, concatenate, check_obs

DATASET = "savas_loi_2018"
BASENAME = "data/{}/unmerged/{}"
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

samples = ["GSM3011853_tils20", "GSM3011854_tils32"]

# merge 10x
adatas = []
for i, sample in enumerate(samples):
    filename = BASENAME.format(DATASET, sample)
    adata = read_10x_mtx(filename, var_names="gene_ids")
    adata.obs["sample"] = sample
    adata.obs["patient"] = str(i)
    adata.obs["origin"] = "tumor_primary"
    adata.obs["replicate"] = 1
    adata.obs["platform"] = "10x_3p_v2"
    adata.obs["tumor_type"] = "BRCA"
    adatas.append(adata)
    
adata = concatenate(adatas, merge_var_cols=["gene_names"])

check_obs(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)

