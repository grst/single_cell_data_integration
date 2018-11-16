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
    script="scripts/02_filter_data/filter_data.Rmd"
  # conda TODO
  output:
    adata=path.join(DATA_PATH_FILTERED, "{dataset}/adata.h5ad"),
    report=path.join(DATA_PATH_FILTERED, "{dataset}/report.html")
  threads: 8
  resources:
    mem_mb=48000
  run:
    param_dict = DATASETS[wildcards.dataset]
    param_dict["doublet_detection"] = True
    param_dict["input_file"] = input.adata
    param_dict["output_file"] = output.adata
    render_rmarkdown(input.script, output.report, ROOT, param_dict)
