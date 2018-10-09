rule preprocess_data:
   """read in all datasets from `data` and bring
   them in a consisten matrix format (either genes x cells
   full matrix or 10x sparse matrix"""
   input:
      "results/book/01_preprocess_data.html"
   output:
      "results/data_processed/datasets.txt",


rule preprocess_data_scanorama:
   """use the scanorama preprocessor to generate `npz` files
   of the data"""
   input:
      "results/data_processed/datasets.txt"
   shell:
      """
      python scanorama/bin/process.py {input}
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
      Rscript -e "rmarkdown::render('{input}', output_file='{output}', output_format=rmarkdown::html_document())"
      """
