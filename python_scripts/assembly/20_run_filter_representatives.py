import os, sys
import pandas as pd
from Bio import SeqIO

def main():

    assembly_filename = os.path.abspath(sys.argv[1])
    assert assembly_filename.endswith(".fasta")
   
    out_filename = os.path.abspath(sys.argv[2])
    assert out_filename.endswith(".fasta")
   
    cluster_reps_filename = os.path.abspath(sys.argv[3])
    assert cluster_reps_filename.endswith(".tsv")

    cluster_reps_df = pd.read_csv(cluster_reps_filename, sep='\t')
    cluster_reps_set = set(cluster_reps_df["cluster"])
	
    to_write = list()
    with open(assembly_filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in cluster_reps_set: to_write.append(record)

    # write records in to_write to .fasta file
    with open(out_filename, "w") as fout:
        SeqIO.write(to_write, fout, "fasta")


if __name__ == "__main__":
    main()

