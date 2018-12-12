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

import scanpy.api as sc
import anndata
import pandas as pd
import os
import sys
sys.path.append("lib")
from scio import concatenate, check_obs, check_var

DATASET = "lambrechts_2018_6149_v2"
OBS_DATA = "pipeline_stages/00_process_fastq/{}/samples.csv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

obs = pd.read_csv(OBS_DATA).drop("batch", axis="columns")
obs = obs.rename({
                        "samples": "samples",
                        "Characteristics[individual]": "patient",
                        "Characteristics[biopsy site]": "origin",
                        "Source Name": "replicate",
                        "Characteristics[sex]": "sex",
                        "Characteristics[age]": "age"
                }, axis="columns")\
         .assign(tumor_type="LUAD", platform="10x_3p_v2")\
        [["samples", "patient", "origin", "replicate", "tumor_type", "platform"]]


dataset_samples = obs["samples"].values


filenames = ["data/{}/{}/raw_gene_bc_matrices_h5.h5".format(DATASET, sample)
             for sample in dataset_samples]

adatas = [sc.read_10x_h5(filename, genome="GRCh38") for filename in filenames]

adatas2 = []
for adata, sample in zip(adatas, dataset_samples):
    # adata.var['gene_symbols'] = adata.var_names
    # adata.var.set_index("gene_ids", inplace=True)
    duplicated = adata.var_names.duplicated()
    print("Removing {} gene symbols because they are duplicated".format(sum(duplicated)))
    adata = adata[:, ~duplicated].copy()
    adata.obs['samples'] = sample
    adatas2.append(adata)

adata = concatenate(adatas2, merge_var_cols=["gene_ids"])
adata.obs = adata.obs.join(obs.set_index("samples"), on="samples", how="left")

adata.obs["dataset"] = DATASET

check_obs(adata)
check_var(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression="lzf")
adata.write_csvs(OUTPUT_DIR)
