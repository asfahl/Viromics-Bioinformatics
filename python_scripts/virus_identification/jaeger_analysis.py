import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.SeqIO.FastaIO import SimpleFastaParser

jaeger_default = pd.read_csv('./1.3_virus_identification/10_jaeger/results_jaeger/assembly/assembly_default_jaeger.tsv', sep='\t', header=0)

with open('./1.2_assembly/10_results_assembly_flye/cross_assembly/assembly.fasta') as fasta_file:  # Will close handle cleanly
    identifiers = []
    lengths = []
    for title, sequence in SimpleFastaParser(fasta_file):
        identifiers.append(title.split(None, 1)[0])  # First word is ID
        lengths.append(sequence)

phages = jaeger_default[jaeger_default["prediction"] == "phage"]

longest = phages["length"].max()
shortest = phages["length"].max()

print(f"The longest phage sequence is {longest}bp long \n")
print(f"The shortest phage sequence is {shortest}bp long")

max_contig_ID = phages["contig_ID"][phages["length"]== longest]

index = identifiers.index(max_contig_ID)
print("The longest contig sequence is: \n")
print(sequence(index))

