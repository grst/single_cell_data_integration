import pysam
import re
import csv
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from collections import defaultdict
import sys
import itertools

#This function fills in a dict with readids
#and their corresponding cell and umi barcodes until it finds
#a specific read id


def parse_barcodes(fastq_parser, query_name, read_barcodes, barcodes_struct):
	for fastq_R1 in fastq_parser:
		read_barcodes[fastq_R1.id]['XC'] = str(fastq_R1.seq)[barcodes_struct['BC_start']:barcodes_struct['BC_end']]
		read_barcodes[fastq_R1.id]['XM'] = str(fastq_R1.seq)[barcodes_struct['UMI_start']:barcodes_struct['UMI_end']]
		if(read_barcodes[fastq_R1.id]['XM']==''):
			sys.SystemExit('UMI empty for read {}.\n The barcode is: {}.\nWhole entry is:{}'.format(fastq_R1.id, fastq_R1.seq,fastq_R1))
		if (fastq_R1.id == query_name):
			return(fastq_parser,read_barcodes)
	return(fastq_parser,read_barcodes)
	
sample_name = 'BT1249'

# cell_barcode_parser = SeqIO.parse(gzip.open(sample_name + '.R2.fastq.gz', "rt"), "fastq")
# umi_barcode_parser = SeqIO.parse(gzip.open(sample_name + '.R3.fastq.gz', "rt"), "fastq")

# outfile = open(sample_name+'_R1.fastq', 'w')

# for cell_barcode_read in cell_barcode_parser:
# 	umi_barcode_read = next(umi_barcode_parser)
	
# 	new_seq = Seq(str(cell_barcode_read.seq) + str(umi_barcode_read.seq[0:10]))
# 	new_qual = cell_barcode_read.letter_annotations['phred_quality'] + umi_barcode_read.letter_annotations['phred_quality'][0:10]
	
# 	cell_barcode_read.letter_annotations={}
# 	cell_barcode_read.seq = new_seq
# 	cell_barcode_read.letter_annotations['phred_quality']=new_qual
# 	SeqIO.write(cell_barcode_read, outfile, 'fastq')

cell_barcode_file = gzip.open(sample_name + '.R2.fastq.gz', "rt")
umi_barcode_file = gzip.open(sample_name + '.R3.fastq.gz', "rt")

outfile = gzip.open(sample_name+'_R1.fastq.gz', 'wb')

line_number = 1
for cell_barcode_line, umi_barcode_line in itertools.zip_longest(cell_barcode_file,umi_barcode_file):
	if(line_number%2==0):
		outfile.write((cell_barcode_line.strip() + umi_barcode_line[0:10] + '\n').encode())
		line_number +=1
		continue
	elif(line_number%4==0):
		outfile.write((cell_barcode_line.strip() + umi_barcode_line[0:10] + '\n').encode())
		line_number +=1
		continue
	else:
		outfile.write(cell_barcode_line.encode())
		line_number +=1
	
