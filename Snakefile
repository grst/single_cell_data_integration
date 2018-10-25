"""
SINGLE CELL DATA INTEGRATION PIPELINE

The purpose of the pipeline is to
* read in single cell datasets ('process_counts') and bring them in to a consistent format
* apply qc, filtering and preprocessing steps to the indvidual datasets ('data_cleaning')
* integrate the datasets by removing batch effects

The rules for the different pipeline stages can be found in the folder `pipeline_stages`.

Each dataset has a unique identifier assigned which is used for folder and file names.
The datasets and their identifiers are defined in `tables/datasets.tsv`
"""

import pandas as pd
import os.path as path

# define the datasets
DATA_PATH = "results/data_processed/"
DATASETS = pd.read_csv("tables/datasets.tsv", sep="\t")["id"].values


include: "pipeline_stages/01_process_counts.smk"


rule render_rmd:
   """render a single rmarkdown document to it's
   corresponding HTML"""
   input:
      "notebooks/{doc}.Rmd"
   output:
      "results/reports/{doc}.html"
   conda:
      "envs/rmarkdown.yml"
   shell:
      """
      cd notebooks && \\
      Rscript -e "reportsdown::render_reports('{input}', \\
        output_file='../{output}', \\
        output_format=reportsdown::html_document2(), \\
        preview=TRUE)"
      """
