import pandas as pd
import sys

# read in the files
blast = sys.argv[1]
taxa = sys.argv[2]
rafah = sys.argv[3]
out = sys.argv[4]

df_rafah = pd.read_csv(rafah, header=0)
print(df_rafah)

df_bins_gtdb = pd.read_csv(taxa, sep='\t' ,header=0)
print(df_bins_gtdb)

df_blast = pd.read_csv(blast, sep='\t', header=None,
                      names=['query', 'subject', 'pid', 'aln_length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
print(df_blast.head(n=10))


# add MAG taxonomy to the blast results
df_blast['bin'] = df_blast['subject'].str.split('.fa').str[0]
print(df_blast.head(n=10))

df_blast_gtdb_merged = df_blast.merge(df_bins_gtdb[['user_genome','classification']], how='left', left_on='bin', right_on='user_genome')
print(df_blast_gtdb_merged.head(n=10))

# merge blast results with rafah
df_blast_gtdb_rafah_merged = df_blast_gtdb_merged.merge(df_rafah,how='left',left_on='query',right_on='Sequence',suffixes=(None,'_rafah'))
df_blast_gtdb_rafah_merged

# select relevant columns
df_blast_gtdb_rafah_merged=df_blast_gtdb_rafah_merged[['query', 'subject', 'pid', 'aln_length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'bin',	'classification','Predicted_Host', 'Winner_Score']]

# rename the columns to meaningful headers
df_blast_gtdb_rafah_merged = df_blast_gtdb_rafah_merged.rename(columns = {'query':'viral_contig', 'subject':'host_bin_contig','bin':'host_bin','classification':'blast_predicted_host','Predicted_Host':'rafah_predicted_host', 'Winner_Score':'rafah_score'})

print(df_blast_gtdb_rafah_merged)

# write out to a csv
df_blast_gtdb_rafah_merged.to_csv(out,index=False)

