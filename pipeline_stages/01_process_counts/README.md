# Preprocess scripts
Bring all datasets in a consistent format. Every dataset
needs to be saved in `results/data_processed` as a `h5ad` file.

## Metadata
Metadata needs to have consistent column names. Columns in **bold** are mandatory.

* **samples**: unique identifier of each sample ('batch')
* **patient**: unique identifier of a patient
* **origin**, origin of the biopsy. controlled vocabs: `tumor_primary`, `normal_adjacent`, `tumor_edge`, `blood_peripheral`, `lymph_node`.
* **replicate**, same patient and origin, but different biopsy
* **platform**, experimental platform. controlled vocabs: `10x_3p_v1`, `10x_3p_v2`, `smartseq2`, `indrop_v2`, `10x_5p`
* **tumor_type**, tissue of origin. use TCGA cancer identifiers, such as `BRCA`.
* **dataset**, the dataset identifier. Used later to visualize batches.


A consistencyt check for the columns and their contents is implemented in
`lib/scio.py:check_obs`.

