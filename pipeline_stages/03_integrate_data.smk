INTEGRATION_TOOLS = ["scanorama", "harmony", "bbknn", "schelker", "ruvg"]

rule integrate_data:
  """
  apply all batch effect removal tool
  on all datasets.
  """
  input:
    expand("results/data_integrated/{tool}/adata.h5ad", tool=INTEGRATION_TOOLS),
    expand("results/data_integrated/{tool}/quality.html", tool=INTEGRATION_TOOLS)
    "results/data_integrated/none/adata.h5ad",
    "results/data_integrated/none/quality.html"


rule _integrate_none:
  """
  integrate data without batch effect removal.
  This is a dependency for all other batch effect removal
  tools as input data and to detect highly variable genes.
  """
  input:
    expand("results/data_filtered/{dataset}/adata.h5ad", dataset=DATASETS),
    script="scripts/03_integrate_data/none.py"
  output:
    adata="results/data_integrated/none/adata.h5ad"
  script:
    input.script


rule _integrate_batch_effect_removal:
  """
  apply a batch effect removal tool
  """
  input:
    expand("results/data_filtered/{dataset}/adata.h5ad", dataset=DATASETS),
    "results/data_integrated/none/adata.h5ad"
    script="scripts/03_integrate_data/{tool}.Rmd"
  output:
    adata="results/data_integrated/{tool}/adata.h5ad",
    report="results/data_integrated/{tool}/report.html"
  threads: 16
  resources:
    mem_mb=128000
  run:
    param_dict = {"output_file": output.adata}
    render_rmarkdown(input.script, output.report, ROOT, param_dict)


rule _integrate_data_report:
  """
  generate a report to assess the quality of the data integration
  """
  input:
    adata="results/data_integrated/{tool}/adata.h5ad",
    script="scripts/03_integrate_data/report_quality.Rmd"
  output:
    report="results/data_integrated/{tool}/quality.html"
  threads: 8
  resources:
    mem_mb=64000
  run:
    param_dict = {"tool": wildcards.tool,
                  "input_data": input.adata}
    render_rmarkdown(input.script, output.report, ROOT, param_dict)

