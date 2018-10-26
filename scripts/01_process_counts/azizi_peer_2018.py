import pandas as pd
import scanpy as sc
from anndata import AnnData
import os.path

from gene_identifiers import map_to_ensembl

DATASET = "azizi_peer_2018"
COUNT_FILE = "data/{}/GSE114725_raw_counts.csv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

raw_counts = pd.read_csv(COUNT_FILE)

obs = raw_counts[["patient", "tissue", "replicate", "cluster", "cellid"]]
obs = obs.assign(cell_name = ["_".join((str(p), str(i))) for p, i in zip(obs.patient, obs.cellid)])
obs = obs.set_index("cell_name")

mat = raw_counts.iloc[:,5:]

var = pd.DataFrame().assign(gene_symbol = mat.columns).set_index("gene_symbol")

adata = AnnData(mat.values, obs, var)
adata = map_to_ensembl(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)


