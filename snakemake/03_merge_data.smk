"""
MERGE DATA PIPELINE STAGE

In this stage, we will merge the datasets into
a single adata object.

In a next step, batch effect corrections can be applied
to the merged object.
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
    script="pipeline_stages/03_merge_data/merge_and_clean.Rmd"
  output:
    adata="results/data_merged/adata.h5ad",
    report="results/data_merged/report.html"
  threads: 32
#  resources:
#    mem_mb=128000
  conda:
    "../envs/merge_data.yml"
  params:
    root_dir=ROOT,
    rmd_params={
      "output_file": "results/data_merged/adata.h5ad",
      "input_files": ",".join(["results/data_processed/{}/adata.h5ad".format(d) for d in DATASETS])
    }
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"

