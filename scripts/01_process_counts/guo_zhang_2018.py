import pandas as pd
import scanpy as sc
from anndata import AnnData
import os.path

DATASET = "guo_zhang_2018"
COUNT_FILE = "data/{}/GSE99254_NSCLC.TCell.S12346.count.txt.gz".format(DATASET)
OBS_FILE = "data/{}/patient_table.tsv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

obs = pd.read_csv(OBS_FILE, sep="\t")
obs.columns = [c.strip() for c in obs.columns]
obs = obs.set_index("UniqueCell_ID")

raw_counts = pd.read_csv(COUNT_FILE, sep="\t")

var = raw_counts[["geneID", "symbol"]].set_index("geneID")

adata = AnnData(raw_counts.iloc[:, 2:].values.transpose(), obs, var)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)
