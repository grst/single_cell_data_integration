import scipy as sp
import os.path

def write_10x(path, matrix, barcodes, genes):
    """Store a count matrix in 10x compatible format.
    
    10x genomics format consists of three files:
     * a matrix-market sparse representation of the data (`matrix.mtx`)
     * a `barcodes.tsv` that contains the cell barcodes
     * a `genes.tsv` that contains the genes 
    Args:
        path: output directory (will contain the three files)
        matrix: matrix (sparse or dense)
        barcodes: pandas dataFrame with the barcodes
        genes: pandas dataFrame with the genes
    """
    assert matrix.shape[0] == barcodes.shape[0], "nrow of matrix does not match length of barcode data frame"
    assert matrix.shape[1] == genes.shape[0], "ncol of matrix does not match length of gene data frame"
    
    sp.mm.write(os.path.join(path, 'matrix.mtx'), matrix)
    barcodes.to_csv(os.path.join(path, 'barcodes.tsv'), sep="\t", header=False, index=False)
    genes.to_csv(os.path.join(path, 'genes.tsv'), sep="\t", header=False, index=False)
    