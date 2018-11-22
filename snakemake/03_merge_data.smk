"""
MERGE DATA PIPELINE STAGE

In this stage, we will merge the datasets into
a single adata object.

Moreover, we will try to reduce the effect of confounders
and normalize the data. All these preprocessing steps
are required for successful application
of batch effect removal tools.
"""

rule merge_data:
  """
  merge all data and regress out confounders.
  """
  input:
    "results/data_merged/adata.h5ad",
    "results/data_merged/report.html"

rule _merge_data:
  """
  merge all data and regress out confounders.
  """
  input:
    expand("results/data_processed/{dataset}/adata.h5ad", dataset=DATASETS),
    script="scripts/03_merge_data/merge_and_clean.Rmd"
  output:
    adata="results/data_merged/adata.h5ad",
    report="results/data_merged/report.html"
  threads: 32
  resources:
    mem_mb=128000
  run:
    param_dict = {"output_data": output.adata}
    render_rmarkdown(input.script, output.report, ROOT, param_dict)

