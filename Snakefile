"""
SINGLE CELL DATA INTEGRATION PIPELINE

The purpose of the pipeline is to
* read in single cell datasets ('process_counts') and bring them in to a consistent format
* apply qc, filtering and preprocessing steps to the indvidual datasets ('data_cleaning')
* integrate the datasets by removing batch effects

The rules for the different pipeline stages can be found in the folder `snakemake`.

Each dataset has a unique identifier assigned which is used for folder and file names.
The datasets and their identifiers are defined in `tables/datasets.tsv`
"""

import pandas as pd
import os.path as path
from lib.snakemaketools import render_rmarkdown

ROOT = path.abspath(path.dirname(workflow.snakefile))

# define the datasets
DATA_PATH = "results/data_processed/"
DATA_PATH_FILTERED = "results/data_filtered"
# DATASETS: dict-like {'zheng_zhang_2017': {'min_genes': 200, 'max_genes':6000, ...}, ...}
DATASETS = pd.read_csv("tables/datasets.tsv", sep="\t", index_col=0).to_dict(orient="index")


include: "snakemake/01_process_counts.smk"
include: "snakemake/02_filter_data.smk"
include: "snakemake/03_merge_data.smk"
# include: "snakemake/04_remove_batch_effects.smk"



