#!/usr/bin/python3

from Bio import SeqIO
import sys 

# read files in as variables
input_gbk = sys.argv[1]
output_gbk = sys.argv[2]
contig = sys.argv[3]

# parse the genbank and write output
for record in SeqIO.parse(input_gbk, "genbank"):
   # print(record.id)
   if record.id == contig:
       SeqIO.write(record, output_gbk, "genbank")

