#!/usr/bin/python3

from Bio import SeqIO
import sys 

# read files in as variables
input_gbk = sys.argv[1]
output_gbk = sys.argv[2]
# set contigs to be extracted
# contigs with few predicted virus hallmark genes per predicted genes
# contigs with highest difference between predicted genes and predicted virus hallmark gene count
contig_12 = ["contig_400", "contig_1459", "contig_1460", "contig_1169", "contig_820", "contig_1062", "contig_200", "contig_58", "contig_267", "contig_1096", "contig_815", "contig_644"]
contig_30 = ["contig_114", "contig_1393", "contig_1160", "contig_103", "contig_185", "contig_1298", "contig_200", "contig_92", "contig_1460", "contig_1552", "contig_1439", "contig_1323", "contig_1433", "contig_1033", "contig_1554", "contig_200", "contig_58", "contig_92", "contig_601", "contig_825", "contig_1096", "contig_644", "contig_1131", "contig_1298", "contig_1323", "contig_1393", "contig_815", "contig_1062", "contig_1170", "contig_267"]


records = []
# parse the genbank and write output
for record in SeqIO.parse(input_gbk, "genbank"):
    # print(record.id)
    if record.id in contig_30:
        records.append(record)
        print(record.id)

print(SeqIO.write(records, output_gbk, "genbank"))

