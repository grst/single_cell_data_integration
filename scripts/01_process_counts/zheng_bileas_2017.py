import pandas as pd
import scanpy.api as sc
from anndata import AnnData
import os.path

DATASET = "zheng_bileas_2017"
COUNT_FILE = "data/{}/".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

adata = sc.read_10x_mtx(COUNT_FILE, var_names="gene_ids") # use ENSG

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)

