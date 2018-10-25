import pandas as pd
import scanpy.api as sc
from anndata import AnnData 
import os.path

DATASET = "savas_loi_2018"
COUNT_FILE = "data/{}/GSE110686_tils20+32_matrix.mtx".format(DATASET)
BARCODES_FILE = "data/{}/GSE110686_tils20+32_barcodes.tsv".format(DATASET)
GENES_FILE = "data/{}/GSE110686_tils20+32_genes.tsv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

adata = sc.read(COUNT_FILE).T

barcodes = pd.read_csv(BARCODES_FILE, sep="\t", header=None)
adata.obs_names = barcodes[0]

genes = pd.read_csv(GENES_FILE, sep="\t", header=None)
adata.var_names = genes[0] #ENSG
adata.var['gene_symbol'] = genes[1].values

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)

