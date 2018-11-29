import pandas as pd
import scanpy as sc
from anndata import AnnData
import os.path
import sys
sys.path.append("lib")
from scio import check_obs, check_var

from gene_identifiers import map_to_ensembl

DATASET = "azizi_peer_2018"
COUNT_FILE = "data/{}/GSE114725_raw_counts.csv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

raw_counts = pd.read_csv(COUNT_FILE)

obs = raw_counts[["patient", "tissue", "replicate", "cluster", "cellid"]]

obs = obs.assign(cell_name = ["_".join((str(p), str(i)))
                             for p, i in zip(obs.patient, obs.cellid)])\
        .assign(origin = "tumor_primary")\
        .assign(platform = "indrop_v2")\
        .assign(tumor_type = "BRCA")

obs = obs.assign(sample = obs[["patient", "origin", "replicate"]].apply(lambda x: "_".join([str(k) for k in x]), axis=1))

obs = obs.set_index("cell_name")

mat = raw_counts.iloc[:,5:]

var = pd.DataFrame().assign(gene_symbol = mat.columns).set_index("gene_symbol")

adata = AnnData(mat.values, obs, var)
adata = map_to_ensembl(adata)

adata.obs["dataset"] = DATASET

check_obs(adata)
check_var(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)

