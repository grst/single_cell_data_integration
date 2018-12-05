from snakemake.io import Namedlist

# tools that operate on the raw gene expression (before cleaning)
INTEGRATION_TOOLS_RAW = [
  "schelker"
]
# tools that operate on the processed, scaled data or on embeddings (e.g. PCA)
INTEGRATION_TOOLS_PROCESSED = ["harmony", "scanorama", "bbknn"]


rule remove_confounders:
  """
  Remove confounders (batch effects, cell cycle, number of genes etc.
  Triggers the ruls for both batch effect removal tools that operate on the
  raw gene expression data (called before scaling the data) and batch
  effect removal tools that are applied as the last processing step
  after scaling and regressing out confounders. """
  input:
    expand("results/data_integrated/final/{tool}/adata.h5ad", tool=INTEGRATION_TOOLS_RAW +
      INTEGRATION_TOOLS_PROCESSED),
    expand("results/data_integrated/final/{tool}/visualize_results.html", tool=INTEGRATION_TOOLS_RAW +
      INTEGRATION_TOOLS_PROCESSED)

rule clean:
  """
  Triggers scaling and regressing out confounders
  as a preprocessing step for batch effect
  removal tools that operate on processed data or
  embeddings.
  """
  input:
    "results/data_integrated/_cleaned/adata.h5ad",
    "results/data_integrated/_cleaned/regress_out_confounders.html"


# This is to work around the issue that snakemake does not
# support input/output in parameters when using a wrapper script.
# In order not to have redundant paths (which is highly probable
# to cause mistakes) I store the paths in the dictionary
# and reference it later.
_integrate_raw_in = {
    "adata" :"results/data_merged/adata.h5ad",
    "notebook" : "pipeline_stages/04_remove_confounders/batch_effect_removal/{tool}.Rmd"
}
_integrate_raw_out = {
    "adata":"results/data_integrated/_integrated/{tool}/adata.h5ad",
    "report":"results/data_integrated/_integrated/{tool}/remove_batch_effects.html"
}
rule _integrate_raw:
  """
  Run batch effect removal tools that operate on the raw gene expression data
  """
  input:
    adata=_integrate_raw_in['adata'],
    notebook=_integrate_raw_in['notebook']
  output:
    adata=_integrate_raw_out['adata'],
    report=_integrate_raw_out['report']
  wildcard_constraints:
    tool="|".join(INTEGRATION_TOOLS_RAW)
  threads: 8
  conda:
    "../envs/{tool}.yml"
  params:
    root_dir=ROOT,
    nb_params = lambda wildcards: dict(
      input_file=_integrate_raw_in['adata'],
      output_file=_integrate_raw_out['adata'].format(tool=wildcards.tool)
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


_clean_raw_in={
    "adata" :"results/data_merged/adata.h5ad",
    "notebook" : "pipeline_stages/04_remove_confounders/regress_out_confounders.Rmd"
}
_clean_raw_out={
    "adata" : "results/data_integrated/_cleaned/adata.h5ad",
    "report" : "results/data_integrated/_cleaned/regress_out_confounders.html"
}
rule _clean_raw:
  """
  Run the 'scaling/regress_out' step
  on the raw gene expression data
  """
  input:
    adata=_clean_raw_in['adata'],
    notebook=_clean_raw_in['notebook']
  output:
    adata=_clean_raw_out['adata'],
    report=_clean_raw_out['report']
  threads: 32
  conda:
    "../envs/merge_data.yml"
  params:
    root_dir=ROOT,
    nb_params = lambda wildcards: dict(
      input_file=_clean_raw_in['adata'],
      output_file=_clean_raw_out['adata']
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


_clean_integrated_in = {
    "adata":"results/data_integrated/_integrated/{tool}/adata.h5ad",
    "notebook":"pipeline_stages/04_remove_confounders/regress_out_confounders.Rmd"
}
_clean_integrated_out = {
    "adata":"results/data_integrated/final/{tool}/adata.h5ad",
    "report":"results/data_integrated/final/{tool}/regress_out_confounders.html"
}
rule _clean_integrated:
  """
  Run the 'scaling/regress_out' step on
  the data integrated by tools that operate on the raw gene
  expression data.
  """
  input:
    adata=_clean_integrated_in['adata'],
    notebook=_clean_integrated_in['notebook']
  output:
    adata=_clean_integrated_out['adata'],
    report=_clean_integrated_out['report']
  wildcard_constraints:
    tool="|".join(INTEGRATION_TOOLS_RAW)
  threads: 16
  conda:
    "../envs/merge_data.yml"
  params:
    root_dir=ROOT,
    nb_params = lambda wildcards: dict(
      input_file=_clean_integrated_in['adata'].format(tool=wildcards.tool),
      output_file=_clean_integrated_out['adata'].format(tool=wildcards.tool)
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


_integrate_cleaned_in = {
    "adata":"results/data_integrated/_cleaned/adata.h5ad",
    "notebook":"pipeline_stages/04_remove_confounders/batch_effect_removal/{tool}.Rmd"
}
_integrate_cleaned_out = {
    "adata":"results/data_integrated/final/{tool}/adata.h5ad",
    "report":"results/data_integrated/final/{tool}/remove_batch_effects.html"
}
rule _integrate_cleaned:
  """
  Run batch effect removal tools
  that operate on the 'cleaned' data.
  """
  input:
    adata=_integrate_cleaned_in['adata'],
    notebook=_integrate_cleaned_in['notebook']
  output:
    adata=_integrate_cleaned_out['adata'],
    report=_integrate_cleaned_out['report']
  wildcard_constraints:
    tool="|".join(INTEGRATION_TOOLS_PROCESSED)
  threads: 8
  conda:
    "../envs/{tool}.yml"
  params:
    root_dir=ROOT,
    nb_params = lambda wildcards: dict(
      input_file=_integrate_cleaned_in['adata'],
      output_file=_integrate_cleaned_out['adata'].format(tool=wildcards.tool)
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


_visualize_results_in = dict(
    adata="results/data_integrated/final/{tool}/adata.h5ad",
    notebook="pipeline_stages/04_remove_confounders/visualize_result.Rmd"
)
_visualize_results_out = {
    "report" : "results/data_integrated/final/{tool}/visualize_results.html"
}
rule _visualize_results:
  """
  Generate a consistent report for each tool
  to assess the quality of integration.
  """
  input:
    adata=_visualize_results_in['adata'],
    notebook=_visualize_results_in['notebook']
  output:
    report=_visualize_results_out['report']
  threads: 8
  conda:
    "../envs/remove_batch_effects.yml"
  params:
    root_dir=ROOT,
    nb_params = lambda wildcards: dict(
      input_file=_visualize_results_in['adata'].format(tool=wildcards.tool),
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


