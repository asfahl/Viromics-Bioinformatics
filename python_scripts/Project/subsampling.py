#!/usr/bin/python3

from Bio import SeqIO
import sys 

# read files in as variables
input_gbk = sys.argv[1]
output_gbk = sys.argv[2]
# set contigs to be extracted
# contigs with few predicted virus hallmark genes per predicted genes
# contigs with highest difference between predicted genes and predicted virus hallmark gene count
contig = ["contig_400", "contig_1459", "contig_1460", "contig_1169", "contig_820", "contig_1062", "contig_200", "contig_58", "contig_267", "contig_1096", "contig_815", "contig_644"]

records = []
# parse the genbank and write output
for record in SeqIO.parse(input_gbk, "genbank"):
    # print(record.id)
    if record.id in contig:
        records.append(record)
        print(record.id)

print(SeqIO.write(records, output_gbk, "genbank"))

