#!/bin/bash
snakemake --use-conda --conda-prefix=$SCRATCH/conda --cluster 'sbatch -t 600 --clusters=serial' --local-cores 4 --jobs 80 --rerun-incomplete concat_files
