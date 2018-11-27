#!/bin/bash
# snakemake --use-conda --conda-prefix /gpfs/scratch/pn69di/ga75tit2/conda  --cluster 'sbatch -t 600 --clusters=serial' --local-cores 4 --jobs 20 --rerun-incomplete download_fastq
snakemake --use-conda --conda-prefix /gpfs/scratch/pn69di/ga75tit2/conda  --cluster 'sbatch -t 600 --clusters=serial' --local-cores 4 --jobs 5 --rerun-incomplete merge_fastq
