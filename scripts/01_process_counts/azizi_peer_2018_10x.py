import pandas as pd
import scanpy as sc
from anndata import AnnData
import os.path
from gene_identifiers import map_to_ensembl

DATASET = "azizi_peer_2018_10x"
COUNT_FILE = "data/{}/GSE114724_raw_counts.tsv".format(DATASET)
OUTPUT_DIR = "results/data_processed/{}/".format(DATASET)

patients = pd.read_csv(COUNT_FILE, sep="\t", nrows=1, header=None)
patients = patients.values[0,1:]
raw_counts = pd.read_csv(COUNT_FILE, sep="\t", skiprows=1)
raw_counts = raw_counts.rename({"Sample_Barcodes": "gene_symbol"}, axis="columns")
raw_counts = raw_counts.set_index("gene_symbol")
barcodes = raw_counts.columns.values

# obs = patients / cells
obs = pd.DataFrame().assign(barcode=barcodes, patient=patients)
obs = obs.set_index("barcode")

# var = genes / transcripts
var = pd.DataFrame().assign(gene_symbol = raw_counts.index.values).set_index("gene_symbol")

adata = AnnData(raw_counts.values.transpose(), obs, var)
adata = map_to_ensembl(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression='lzf')
adata.write_csvs(OUTPUT_DIR)
