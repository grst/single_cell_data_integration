"""
Data filtering pipeline stage.

Filter every dataset individually before merging based on QC criteria.
These include
* minimum and maximum number of detected genes
* fraction of mitochondrial reads
* doublet detection.

In this pipeline stage, the data *as kept as raw counts*
with *no normalization* as this is required for successful merging.
"""

rule filter_data:
  """
  run the cleaning/filtering script on all dataset
  with the parameters specified in the
  datasets.tsv table
  """
  input:
    "tables/datasets.tsv",
    expand(path.join(DATA_PATH_FILTERED, "{dataset}/adata.h5ad"), dataset=DATASETS)


rule _filter_data:
  """
  apply cleaning/filtering script to a dataset
  """
  input:
    adata=path.join(DATA_PATH, "{dataset}/adata.h5ad"),
    script="pipeline_stages/02_filter_data/filter_data.Rmd"
  output:
    adata=path.join(DATA_PATH_FILTERED, "{dataset}/adata.h5ad"),
    report=path.join(DATA_PATH_FILTERED, "{dataset}/report.html")
  conda: "../envs/filter_data.yml"
  threads: 8
  resources:
    mem_mb=48000
  params:
    root_dir=ROOT,
    rmd_params=lambda wildcards: dict(DATASETS[wildcards.dataset],
        input_file=path.join(DATA_PATH, "{}/adata.h5ad".format(wildcards.dataset)),
        output_file=path.join(DATA_PATH_FILTERED, "{}/adata.h5ad".format(wildcards.dataset)),
        doublet_detection=True)
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"
