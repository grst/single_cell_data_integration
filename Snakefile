import pandas as pd
DATA_PATH = "results/data_processed/"
DATASETS = pd.read_csv("tables/datasets.tsv", sep="\t")
DATA_FILES_PROCESSED = []
DATA_FILES_SCANORAMA = []
for i, path, type in DATASETS.itertuples():
  if type == 'tsv':
    DATA_FILES_PROCESSED.append(path + '.txt')
    DATA_FILES_SCANORAMA.append(path + '.npz')
  elif type == "mtx":
    DATA_FILES_PROCESSED.extend([path + '/genes.tsv', path + '/barcodes.tsv', path + '/matrix.mtx'])
    DATA_FILES_SCANORAMA.extend([path + '/tab.genes.txt', path + '/tab.npz'])


rule preprocess_data:
   """read in all datasets from `data` and bring
   them in a consisten matrix format (either genes x cells
   full matrix or 10x sparse matrix"""
   input:
      "results/book/01_preprocess_data.html",
      "tables/datasets.tsv"
   output:
      expand(DATA_PATH + "{data_file}", data_file=DATA_FILES_PROCESSED)


rule scanorama_input_file:
   """generate text file that list the dataset in a scanorama-compatible way. """
   input:
      "tables/datasets.tsv"
   output:
      DATA_PATH + "scanorama_datasets.txt"
   run:
      with open(output[0], 'wt') as f:
        for p in DATASETS.path:
          f.write(DATA_PATH + p + '\n')


rule preprocess_data_scanorama:
   """use the scanorama preprocessor to generate `npz` files
   of the data"""
   input:
      expand(DATA_PATH + "{data_file}", data_file=DATA_FILES_PROCESSED),
      scanorama_datasets=DATA_PATH + "scanorama_datasets.txt"
   output:
      expand(DATA_PATH + "{data_file_npz}", data_file_npz=DATA_FILES_SCANORAMA)
   conda:
      "envs/scanorama.yml"
   shell:
      """
      conda develop scanorama
      python scanorama/bin/process.py {input.scanorama_datasets}
      """

rule scanorama:
   """run scanorama to integrate the single cell datasets. """
   input:
      expand(DATA_PATH + "{data_file_npz}", data_file_npz=DATA_FILES_SCANORAMA),
      scanorama_datasets=DATA_PATH + "scanorama_datasets.txt"
   output:
      "results/scanorama/datasets_dimred.npz",
      "results/scanorama/datasets.npz",
      "results/scanorama/genes.npz",
      "results/scanorama/tsne.npz"
   conda:
      "envs/scanorama.yml"
   shell:
      """
      conda develop scanorama
      python scanorama/bin/run_panorama.py {input.scanorama_datasets} results/scanorama
      """



rule render_rmd:
   """render a single rmarkdown document to it's
   corresponding HTML"""
   input:
      "notebooks/{doc}.Rmd"
   output:
      "results/book/{doc}.html"
   conda:
      "envs/rmarkdown.yml"
   shell:
      """
      cd notebooks && \\
      Rscript -e "bookdown::render_book('{input}', \\
        output_file='../{output}', \\
        output_format=bookdown::html_document2(), \\
        preview=TRUE)"
      """
