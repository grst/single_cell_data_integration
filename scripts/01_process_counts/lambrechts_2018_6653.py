import scanpy.api as sc
import anndata
import pandas as pd
import os
import sys
sys.path.append("lib")
from scio import concatenate

DATASET = "lambrechts_2018_6653"
OBS_DATA = "raw/{}/samples.csv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

obs = pd.read_csv(OBS_DATA).drop("batch", axis="columns")
obs = obs.rename({
                        "samples": "sample",
                        "Characteristics[individual]": "patient",
                        "Characteristics[sampling site]": "origin",
                        "Source Name": "replicate",
                        "Characteristics[sex]": "sex",
                        "Characteristics[age]": "age"
                }, axis="columns")\
         .assign(tumor_type="LUAD", platform="10x_3p_v2")\
        [["sample", "patient", "origin", "replicate", "tumor_type", "platform"]]

origin_map = {
    "tumour edge": "tumor_edge",
    "tumour middle (in between core and edge sample)": "tumor_primary",
    "tumour core": "tumor_primary",
    "normal tissue adjacent to tumour": "normal_adjacent"
}

obs = obs.assign(origin = obs["origin"].apply(lambda x: origin_map[x]))\
         .assign(replicate = obs["replicate"].apply(lambda x: x[-1]))


dataset_samples = obs["sample"].values


filenames = ["raw/{}/data/cellranger/{}/outs/raw_gene_bc_matrices_h5.h5".format(DATASET, sample) 
             for sample in dataset_samples]

adatas = [sc.read_10x_h5(filename, genome="GRCh38") for filename in filenames]

for adata, sample in zip(adatas, dataset_samples):
    adata.var['gene_name'] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)
    adata.obs['sample'] = sample
    
adata = concatenate(adatas, merge_var_cols=["gene_name"])
adata.obs = adata.obs.join(obs.set_index("sample"), on="sample", how="left")

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression="lzf")
adata.write_csvs(OUTPUT_DIR)