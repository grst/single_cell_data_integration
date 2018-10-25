# ---
# jupyter:
#   jupytext_format_version: '1.2'
#   kernelspec:
#     display_name: Python [conda env:single_cell_integration]
#     language: python
#     name: conda-env-single_cell_integration-py
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.6.6
# ---

import pandas as pd
import scanpy as sc
from anndata import AnnData
import os.path

DATASET = "azizi_peer_2018"
COUNT_FILE = "../../data/{}/GSE114725_raw_counts.csv".format(DATASET)
OUTPUT_DIR = "../../results/data_processed/{}/".format(DATASET)

raw_counts = pd.read_csv(COUNT_FILE)

obs = raw_counts[["patient", "tissue", "replicate", "cluster", "cellid"]]
obs = meta.assign(cell_name = ["_".join((str(p), str(i))) for p, i in zip(meta.patient, meta.cellid)])
obs = meta.set_index("cell_name")

mat = raw_counts.iloc[:,5:]

var = pd.DataFrame().assign(gene_symbol = mat.columns).set_index("gene_symbol")

adata = AnnData(mat.values, obs, var)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)


