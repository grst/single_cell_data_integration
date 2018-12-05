#!/usr/bin/env python3
"""preprocess 10x' chromium v1 chemistry read files for the use
with dropSeqPipeself.

Usage:
  preprocess_10x_3p_v1.py I1_FILE RA_FILE R1_OUTPUT R2_OUTPUT
"""


from Bio.SeqIO.QualityIO import FastqGeneralIterator
from docopt import docopt
from os.path import abspath, expanduser
import os.path as path
from glob import glob
from shutil import copyfile
import gzip


def main(i1_file, ra_file, r1_output, r2_output):
    with gzip.open(i1_file, 'rt') as i1_h:
        with gzip.open(ra_file, 'rt') as ra_h:
            with gzip.open(r1_output, 'wt') as r1_h:
                with gzip.open(r2_output, 'wt') as r2_h:
                    i1_it = FastqGeneralIterator(i1_h)
                    ra_it = FastqGeneralIterator(ra_h)
                    for barcode_record in i1_it:
                        read_record = next(ra_it)
                        umi_record = next(ra_it)
                        assert barcode_record[0].split()[0] == \
                            read_record[0].split()[0] == \
                            umi_record[0].split()[0], "record headers match"
                        # write read
                        r2_h.write("@" + read_record[0] + "\n")
                        r2_h.write(read_record[1] + "\n")
                        r2_h.write("+\n")
                        r2_h.write(read_record[2] + "\n")
                        # write barcode + umi
                        r1_h.write("@" + read_record[0] + "\n")
                        r1_h.write(barcode_record[1] + umi_record[1] + "\n")
                        r1_h.write("+\n")
                        r1_h.write(barcode_record[2] + umi_record[2] + "\n")

if __name__ == '__main__':
    arguments = docopt(__doc__)
    main(arguments['I1_FILE'], arguments['RA_FILE'], arguments['R1_OUTPUT'], arguments["R2_OUTPUT"])
