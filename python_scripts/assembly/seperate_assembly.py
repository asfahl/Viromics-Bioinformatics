import os, sys
from Bio import SeqIO

def main():

   # access parameters passed to your script
   assembly = os.path.abspath(sys.argv[1])

   # throw an error if the statement behind assert is not true
   assert assembly.endswith(".fasta")
   
   out_fasta = os.path.abspath(sys.argv[2])
   assert out_fasta.endswith(".fasta")
   
   sample_id = sys.argv[3]
   
   # modify names of the scaffols and store in to_write list
   to_write = list()
   with open(assembly) as handle:
       for record in SeqIO.parse(handle, "fasta"):
           # make sure not to do this twice, if run several times. Its not necessary.
           if record.id.startswith(sample_id): continue

           # adjust the id of the record object
           record.id = f"{sample_id}_{record.id}"

           # append the object to the list for writing back to a file later.
           to_write.append(record)

   # write records in to_write to .fasta file
   with open(out_fasta, "w") as fout:
       SeqIO.write(to_write, fout, "fasta")


if __name__ == "__main__":
   main()

