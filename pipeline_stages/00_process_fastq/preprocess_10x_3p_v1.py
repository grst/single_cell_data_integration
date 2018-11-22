#!/usr/bin/env python3
"""preprocess 10x' chromium v1 chemistry read files for the use
with dropSeqPipeself.

Usage:
  preprocess_10x_3p_v1.py BASENAME INPUT_DIR OUTPUT_DIR

Arguments:
  BASENAME      part of the filename that uniquely identifies a sampleself.
                e.g. `sample_1`, if the read files are
                `sample_1.{R1,R2,R3}.fastq`
  INPUT_DIR     directory that are searched for the fastq files
  OUTPUT_DIR    output directory, where the processed files will be saved.
"""


from Bio.SeqIO.QualityIO import FastqGeneralIterator
from docopt import docopt
from os.path import abspath, expanduser
import os.path as path
from glob import glob
from shutil import copyfile
import gzip

CONFIG = {
    'umi-length': 10,
    'barcode-length': 14
}


def main(basename, input_dir, output_dir):
    input_dir = abspath(expanduser(input_dir))
    output_dir = abspath(expanduser(output_dir))
    output_r1 = path.join(output_dir, basename + "_R1.fastq.gz")
    output_r2 = path.join(output_dir, basename + "_R2.fastq.gz")
    r1_file = glob(input_dir + '/*{}*R1*.fastq.gz'.format(basename))
    r2_file = glob(input_dir + '/*{}*R2*.fastq.gz'.format(basename))
    r3_file = glob(input_dir + '/*{}*R3*.fastq.gz'.format(basename))

    if(len(r1_file) != 1):
        raise Exception("More than 1 R1 file found")
    if(len(r2_file) != 1):
        raise Exception("More than 1 R2 file found")
    if(len(r3_file) != 1):
        raise Exception("More than 1 R3 file found")

    with gzip.open(r2_file[0], 'rt') as barcode_handle:
        with gzip.open(r3_file[0], 'rt') as umi_handle:
            with gzip.open(output_r1, 'wt') as out_handle:
                for barcode_record, umi_record in zip(
                                        FastqGeneralIterator(barcode_handle),
                                        FastqGeneralIterator(umi_handle)):
                    assert barcode_record[0].split()[0] == \
                        umi_record[0].split()[0], "record titles match"
                    seq_string = barcode_record[1][:CONFIG['barcode-length']] + umi_record[1][:CONFIG['umi-length']]
                    qual_string = barcode_record[2][:CONFIG['barcode-length']] + umi_record[2][:CONFIG['umi-length']]
                    out_handle.write("@" + barcode_record[0] + "\n")
                    out_handle.write(seq_string + "\n")
                    out_handle.write("+\n")
                    out_handle.write(qual_string + "\n")

    copyfile(r1_file[0], output_r2)


if __name__ == '__main__':
    arguments = docopt(__doc__)
    main(arguments['BASENAME'], arguments['INPUT_DIR'], arguments['OUTPUT_DIR'])
