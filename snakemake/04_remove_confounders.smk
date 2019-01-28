"""
Remove
  * coufounding factors, such as fraction of mito-genes or
    reads per cell
  * batch effects, i.e. the influence of datasets using
    different tools
  * normalize and scale the data -> prepare for PCA/UMAP/downstream
    analysis.

"""


from snakemake.io import Namedlist

# tools that operate on the raw gene expression (before cleaning)
INTEGRATION_TOOLS_RAW = [
  "schelker"
]
# tools that operate on the processed, scaled data or on embeddings (e.g. PCA)
INTEGRATION_TOOLS_PROCESSED = ["harmony", "scanorama", "bbknn", "harmony-patient", "bbknn-patient"]


rule remove_confounders:
  """
  Remove confounders (batch effects, cell cycle, number of genes etc.
  Triggers the ruls for both batch effect removal tools that operate on the
  raw gene expression data (called before scaling the data) and batch
  effect removal tools that are applied as the last processing step
  after scaling and regressing out confounders. """
  input:
    expand("results/data_integrated/batch_effects_removed/{tool}/adata.h5ad", tool=INTEGRATION_TOOLS_RAW +
      INTEGRATION_TOOLS_PROCESSED),
    expand("results/data_integrated/batch_effects_removed/{tool}/visualize_results.html", tool=INTEGRATION_TOOLS_RAW +
      INTEGRATION_TOOLS_PROCESSED),
    "results/data_integrated/compare_tools/comparison.html"

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

rule annotate_cell_types:
  """
  run the automated cell type annotation on the cleaned
  data.
  """
  input:
    "results/data_integrated/cell_types/cell_types.tsv"


rule _compare_tools:
  """
  Execute the notebook that compares batch effect removal tools
  """
  input:
    expand("results/data_integrated/batch_effects_removed/{tool}/adata.h5ad",
              tool=INTEGRATION_TOOLS_RAW+INTEGRATION_TOOLS_PROCESSED),
    "results/data_integrated/cell_types/cell_types.tsv",
    notebook="pipeline_stages/04_remove_confounders/compare_tools.Rmd"
  output:
    report="results/data_integrated/compare_tools/comparison.html"
  conda:
    "../envs/compare_batch_effect_removal_tools.yml"
  threads: 16
  params:
    root_dir=ROOT,
    nb_params = dict()
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


# This is to work around the issue that snakemake does not
# support input/output in parameters when using a wrapper script.
# In order not to have redundant paths (which is highly probable
# to cause mistakes) I store the paths in the dictionary
# and reference it later.

# Run integration tools that use raw counts (before regressing
# out confounders and scaling).
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
    notebook=_integrate_raw_in['notebook'],
  output:
    adata=_integrate_raw_out['adata'],
    report=_integrate_raw_out['report']
  wildcard_constraints:
    tool="|".join(INTEGRATION_TOOLS_RAW)
  threads: 32   # needed to switch to float64 matrix, needs a lot of memory now.
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


# Regress out confounders and scale/normalize raw data.
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
  on the raw gene expression data. This
  is done before calling harmony/bbknn/scanorama or similar tools.
  """
  input:
    adata=_clean_raw_in['adata'],
    notebook=_clean_raw_in['notebook']
  output:
    adata=_clean_raw_out['adata'],
    report=_clean_raw_out['report']
  threads: 32
  conda:
    "../envs/remove_batch_effects.yml"
  params:
    root_dir=ROOT,
    nb_params = lambda wildcards: dict(
      input_file=_clean_raw_in['adata'],
      output_file=_clean_raw_out['adata']
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


# Annotate cell types
_annotate_ct_in={
    "adata_scaled" : "results/data_integrated/_cleaned/adata.h5ad",
    "adata_raw" : "results/data_merged/adata.h5ad",
    "notebook" : "pipeline_stages/04_remove_confounders/annotate_cell_types.Rmd"
}
_annotate_ct_out={
    "table" : "results/data_integrated/cell_types/cell_types.tsv",
    "report" : "results/data_integrated/cell_types/annotate_cell_types.html"
}
rule _annotate_cell_types:
  """
  Annotate cell types automatically.
  Runs on the cleaned/scaled data, but before
  batch effect correction. (Need this information
  to assess batch effect correction tools).
  """
  input:
    adata_raw=_annotate_ct_in['adata_raw'],
    adata_scaled=_annotate_ct_in['adata_scaled'],
    notebook=_annotate_ct_in['notebook']
  output:
    table=_annotate_ct_out['table'],
    report=_annotate_ct_out['report']
  threads: 16
  conda:
    "../envs/annotate_cell_types.yml"
  params:
    root_dir=ROOT,
    nb_params = lambda wildcards: dict(
      adata_raw=_annotate_ct_in['adata_raw'],
      adata_scaled=_annotate_ct_in['adata_scaled'],
      output_file=_annotate_ct_out['table']
    )
  wrapper:
    "file:snakemake/wrappers/render_rmarkdown"


# Run the 'regress_out_coufounders' step on
# already integrated data (e.g. Schelker's approach)
_clean_integrated_in = {
    "adata":"results/data_integrated/_integrated/{tool}/adata.h5ad",
    "notebook":"pipeline_stages/04_remove_confounders/regress_out_confounders.Rmd"
}
_clean_integrated_out = {
    "adata":"results/data_integrated/batch_effects_removed/{tool}/adata.h5ad",
    "report":"results/data_integrated/batch_effects_removed/{tool}/regress_out_confounders.html"
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


# Run harmony/scanorama/bbknn on the
# already scaled data.
_integrate_cleaned_in = {
    "adata":"results/data_integrated/_cleaned/adata.h5ad",
    "notebook":"pipeline_stages/04_remove_confounders/batch_effect_removal/{tool}.Rmd"
}
_integrate_cleaned_out = {
    "adata":"results/data_integrated/batch_effects_removed/{tool}/adata.h5ad",
    "report":"results/data_integrated/batch_effects_removed/{tool}/remove_batch_effects.html"
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


# Visualize the results of the different integration tools.
_visualize_results_in = dict(
    adata="results/data_integrated/batch_effects_removed/{tool}/adata.h5ad",
    notebook="pipeline_stages/04_remove_confounders/visualize_result.Rmd"
)
_visualize_results_out = {
    "report" : "results/data_integrated/batch_effects_removed/{tool}/visualize_results.html"
}
rule _visualize_results:
  """
  Generate a consistent report for each tool
  to assess the quality of integration.
  """
  input:
    adata=_visualize_results_in['adata'],
    notebook=_visualize_results_in['notebook'],
    cell_types="results/data_integrated/cell_types/cell_types.tsv"
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


