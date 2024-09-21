# This script takes 3 arguments. Run it:
# python your_script_file_name.py path/to/assembly1.fasta path/to/assembly2.fasta box_and_violin.png

import os, sys
from Bio import SeqIO
import matplotlib.pyplot as plt

def main():
    # take the first argument to the script as the filename to the bacterial assembly 
    assembly1_filename = os.path.abspath(sys.argv[1])
    assert assembly1_filename.endswith(".fasta")

    # take the second argument to the script as the filename to the viral assembly
    assembly2_filename = os.path.abspath(sys.argv[2])
    assert assembly2_filename.endswith(".fasta")

    # set the filename of the output PNG file
    out_filename = os.path.abspath(sys.argv[3])
    assert out_filename.endswith(".png")

    # open the assembly files and get the lengths of each contig within one list 
    with open(assembly1_filename) as handle:
        lengths_assembly1 = [len(record.seq) for record in SeqIO.parse(handle, "fasta")]

    with open(assembly2_filename) as handle:
        lengths_assembly2 = [len(record.seq) for record in SeqIO.parse(handle, "fasta")]

    # use pyplot to create a multi panel plot with 2 columns and 1 row. figsize is in inches...
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))

    # make violin plots for both lists of contig lengths
    axs[0].violinplot([lengths_assembly1, lengths_assembly2], showextrema=False)

    # add axis description and title, limit the y range for comparability
    axs[0].set_xticks([1,2], ["phages", "bacteria"])
    axs[0].set_ylim([0,100000])
    axs[0].set_title('Violin plot')

    # make box plots for both lists of contig lengths and add lables
    axs[1].boxplot([lengths_assembly1, lengths_assembly2], showfliers=False)
    axs[1].set_xticks([1,2], ["phages", "bacteria"])
    axs[1].set_title('Box plot')
    
    fig.savefig(out_filename, dpi=200)
    

if __name__ == "__main__":
    main()

