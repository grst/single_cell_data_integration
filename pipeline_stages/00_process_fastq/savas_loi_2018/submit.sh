# savas_loi_2018 was obtained from zenodo as BCL files specified in their manuscript.
# convert to fastq using cellranger 2.2.0 and bcl2fastq v2.20.0.422
# --use-bases-mask was derived from information on https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/bcl2fastq-direct/ and from the RunInfo.xml files.
conda install -c freenome bcl2fastq
/storage/apps/pipelines/cellranger_10x/cellranger-2.2.0/cellranger mkfastq --localcores=32 --localmem=150 --run AGRF_CAGRF14386_CAUC0ANXX  --simple-csv AGRF_CAGRF14386_CAUC0ANXX/simple_SampleSheet.csv --use-bases-mask="Y26n*,I8n,Y101"
/storage/apps/pipelines/cellranger_10x/cellranger-2.2.0/cellranger mkfastq --localcores=32 --localmem=150 --run AGRF_CAGRF14386_CB12LANXX  --simple-csv AGRF_CAGRF14386_CB12LANXX/simple_SampleSheet.csv --use-bases-mask="Y26n*,I8n,Y101"

