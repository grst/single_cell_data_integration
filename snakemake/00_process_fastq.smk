from lib.scio import check_samples_csv

rule check_samples_csv:
  input:
    "tables/datasets.tsv"
  run:
    for ds in DATASETS:
      print("##### Checking samples.csv for {}".format(ds))
      try:
        check_samples_csv("pipeline_stages/00_process_fastq/{}/samples.csv".format(ds))
      except AssertionError as e:
        print(e)
      print("\n")

