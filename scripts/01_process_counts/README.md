# Preprocess scripts
Bring all datasets in a consistent format. Every dataset 
needs to be saved in `results/data_processed` as a `h5ad` file. 

## Metadata
Metadata needs to have consistent column names. Columns in **bold** are mandatory.

* **sample**: unique identifier of each sample ('batch')
* **patient**: unique identifier of a patient
* **origin**, origin of the biopsy. controlled vocabs: `tumor_primary`, `normal_adjacent`, `tumor_edge`
* **replicate**, same patient and origin, but different biopsy
* **platform**, experimental platform. controlled vocabs: `10x_3p`, `10x_3p_v2`, `smartseq2`, `indrop_v2`
* **tumor_type**, tissue of origin. use TCGA cancer identifiers, such as `BRCA`.

