import os, sys
from Bio import SeqIO  # for handling biological sequence data (like DNA, RNA etc.) in Python.

def main():
    # define file paths from the arguments
    assembly_path = os.path.abspath(sys.argv[1])
    assert assembly_path.endswith(".fasta")

    # set an output directory and create it if it does not exist
    out_dir = os.path.abspath(sys.argv[2])
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    # loop through the records in the combined assembly
    with open(assembly_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # set a filename per record
            out_fasta = os.path.join(out_dir, f"{record.id}.fasta")
            # write the record to the file
            with open(out_fasta, "w") as fout:
                SeqIO.write([record], fout, "fasta")
                
if __name__ == "__main__":
    main()