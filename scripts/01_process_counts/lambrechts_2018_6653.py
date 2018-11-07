import scanpy.api as sc
import anndata
import pandas as pd
import os

DATASET = "lambrechts_2018_6653"
OBS_DATA = "raw/{}/samples.csv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

obs = pd.read_csv(OBS_DATA).drop("batch", axis="columns")
dataset_samples = obs.samples.values


filenames = ["raw/{}/data/cellranger/{}/outs/raw_gene_bc_matrices_h5.h5".format(DATASET, sample) 
             for sample in dataset_samples]

adatas = [sc.read_10x_h5(filename, genome="GRCh38") for filename in filenames]

for adata, sample in zip(adatas, dataset_samples):
    adata.var['gene_name'] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)
    adata.obs['sample'] = sample
    
adata = adatas[0].concatenate(adatas[1:])
adata.obs = adata.obs.join(obs.set_index("samples"), on="sample", how="left")

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression="lzf")
adata.write_csvs(OUTPUT_DIR)