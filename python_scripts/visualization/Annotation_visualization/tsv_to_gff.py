import pandas as pd
import gffutils.gffwriter as gw
import gffutils.feature as gf
import os, sys

def main():

	in_file = os.path.abspath(sys.argv[1])
	assert in_file.endswith(".tsv")

	out_file = os.path.abspath(sys.argv[2])
	writer = gw.GFFWriter(out_file)
	
	contig_id = sys.argv[3]
	
	df = pd.read_csv(in_file, sep="\t")
	df = df.loc[df["contig_id"] == contig_id]
	
	for index, row in df.iterrows():
		if row["strand"] == 1: strand = "+"
		else: strand = "-"
		f = gf.Feature(seqid=row["contig_id"], featuretype="CDS", start=row["start"], end=row["end"], strand=strand)
		writer.write_rec(f)
	
	writer.close()

if __name__ == "__main__":
    main()


    