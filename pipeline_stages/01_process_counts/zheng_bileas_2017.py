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
import sys
sys.path.append("lib")
from scio import check_obs, check_var

DATASET = "zheng_bileas_2017"
COUNT_FILE = "data/{}/".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

adata = sc.read_10x_mtx(COUNT_FILE, var_names="gene_symbols")

adata.obs = adata.obs.assign(samples = "1")\
                     .assign(patient = "1")\
                     .assign(origin = "blood_peripheral")\
                     .assign(replicate = "1")\
                     .assign(platform = "10x_3p_v2")\
                     .assign(tumor_type = "PBMC")

adata.obs["dataset"] = DATASET

check_obs(adata)
check_var(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
