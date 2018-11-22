"""
RULES FOR THE 'PROCESS COUNTS' PIPELINE STAGE

At this stage, we will
* read in the single cell datasets
* map identifiers
* generate scanpy AnnData objects
* and store them in the 'processed_data' results folder.


For each dataset, a preprocess script named [UNIQUE_IDENTIFIER].py
needs to be placed in `scripts/01_process_counts/`
"""

rule process_counts:
    """
    run the preprocess scripts on all datasets
    """
    input:
        "tables/datasets.tsv",
        expand(path.join(DATA_PATH, "{dataset}/adata.h5ad"), dataset=DATASETS)


rule _process_counts:
    """
    apply process_counts on a dataset.
    """
    input:
        "scripts/01_process_counts/{dataset}.py"
    conda:
        "../envs/process_counts.yml"
    output:
        path.join(DATA_PATH, "{dataset}/adata.h5ad")
    shell:
        "PYTHONPATH=lib python {input}"

