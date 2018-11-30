INTEGRATION_TOOLS = ["schelker", "scanorama", "bbknn", "harmony"]

rule remove_batch_effects:
  """
  apply all batch effect removal tool
  on all datasets.
  """
  input:
    expand("results/data_integrated/{tool}/adata.h5ad", tool=INTEGRATION_TOOLS),
    expand("results/data_integrated/{tool}/quality.html", tool=INTEGRATION_TOOLS)


rule _integrate_batch_effect_removal:
  """
  apply a batch effect removal tool
  """
  input:
    "results/data_merged/adata.h5ad",
    script="pipeline_stages/04_remove_batch_effects/{tool}.Rmd"
  output:
    adata="results/data_integrated/{tool}/adata.h5ad",
    report="results/data_integrated/{tool}/report.html"
  threads: 16
  resources:
    mem_mb=128000
  conda:
    "../envs/{tool}.yml"
  params:
    root_dir=ROOT,
    rmd_params = lambda wildcards: dict(
      output_file="results/data_integrated/{}/adata.h5ad".format(wildcards.tool),
      input_file="results/data_merged/adata.h5ad"
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


rule _remove_batch_effects_report:
  """
  generate a report to assess the quality of the data integration
  """
  input:
    adata="results/data_integrated/{tool}/adata.h5ad",
    script="pipeline_stages/04_remove_batch_effects/visualize.Rmd"
  output:
    report="results/data_integrated/{tool}/quality.html"
  threads: 8
  resources:
    mem_mb=64000
  conda:
    "../envs/remove_batch_effects.yml"
  params:
    root_dir=ROOT,
    rmd_params = lambda wildcards: dict(
      input_file="results/data_integrated/{}/adata.h5ad".format(wildcards.tool)
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"
