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
    adata = read_10x_mtx(filename, var_names="gene_ids")
    adata.obs["sample"] = row["sample"]
    adatas.append(adata)

adata = concatenate(adatas, merge_var_cols=["gene_symbols"])

adata.obs = adata.obs.join(obs.set_index('sample'), on="sample", how="left")

adata.obs["dataset"] = DATASET

check_obs(adata)
check_var(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression="lzf")
adata.write_csvs(OUTPUT_DIR)
