import os, sys
import pandas as pd
from Bio import SeqIO

def main():

    assembly_path = os.path.abspath(sys.argv[1])
    assert assembly_path.endswith(".fasta")

    jaeger_results_path = os.path.abspath(sys.argv[2])
    assert jaeger_results_path.endswith(".tsv")
    
    checkv_results_path = os.path.abspath(sys.argv[3])
    assert checkv_results_path.endswith(".tsv")
    
    out_fasta = os.path.abspath(sys.argv[4])
    assert out_fasta.endswith(".fasta")

    # read the tsv files as pandas dataframes
    jaeger_df = pd.read_csv(jaeger_results_path, sep='\t')
    checkv_df = pd.read_csv(checkv_results_path, sep='\t')
	
    # collect the sets of contigs which stick to our selection cutoffs
    jaeger_selection = {row['contig_id'] for index, row in jaeger_df.iterrows() if row['prediction'] == 'phage'}
    checkv_selection = {row['contig_id'] for index, row in checkv_df.iterrows() if row['completeness'] > 50 and row['contamination'] < 5}
	
    # use set operation union to get the contigs in the jaeger set AND in the checkv set
    joint_selection = jaeger_selection.intersection(checkv_selection)
	
    # print some numbers
    print(f"Number of input contigs: {len(jaeger_df.index)}, selected by jaeger: {len(jaeger_selection)}, selected by checkv: {len(checkv_selection)}, joint selection: {len(joint_selection)}")
	
    # define list of records to keep and fill it by comparing the contig id of each record to the joint set of selected contigs
    out_records = []
    with open(assembly_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in joint_selection: out_records.append(record)
    
    # write the selected records into a new file
    with open(out_fasta, "w") as fout:
        SeqIO.write(out_records, fout, "fasta")

if __name__ == "__main__":
    main()

