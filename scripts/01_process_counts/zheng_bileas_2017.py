import pandas as pd
import scanpy.api as sc
from anndata import AnnData
import os.path
import sys
sys.path.append("lib")
from scio import check_obs

DATASET = "zheng_bileas_2017"
COUNT_FILE = "data/{}/".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

adata = sc.read_10x_mtx(COUNT_FILE, var_names="gene_ids") # use ENSG

adata.obs = adata.obs.assign(sample = "1")\
                     .assign(patient = "1")\
                     .assign(origin = "blood_peripheral")\
                     .assign(replicate = "1")\
                     .assign(platform = "10x_3p_v2")\
                     .assign(tumor_type = "PBMC")

check_obs(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)

