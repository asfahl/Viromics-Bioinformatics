# %%
import pandas as pd
import functools as ft

assembly_filename = "/home/qe76qox/Viromics-Bioinformatics/1.2_assembly/10_results_assembly_flye/cross_assembly/assembly_info.txt"
jaeger_filename = "/home/qe76qox/Viromics-Bioinformatics/1.3_virus_identification/10_jaeger/results_jaeger/assembly/assembly_default_jaeger.tsv"
checkv_filename = "/home/qe76qox/Viromics-Bioinformatics/1.3_virus_identification/20_results_assessment_checkv/cross_assembly/quality_summary.tsv"
pharokka_filename = "/home/qe76qox/Viromics-Bioinformatics/1.4_annotation/10_pharokka/pharokka_cds_functions.tsv"
genomad_filename = "/home/qe76qox/Viromics-Bioinformatics/2.2_taxonomy/20_genomad_results/assembly_annotate/assembly_genes.tsv"
genomad_taxonomy_filename = "/home/qe76qox/Viromics-Bioinformatics/2.2_taxonomy/20_genomad_results/assembly_annotate/assembly_taxonomy.tsv"
vcontact_filename = "/home/qe76qox/Viromics-Bioinformatics/2.2_taxonomy/10_vcontact_results/final_assignments.csv"
out_filename = "/home/qe76qox/Viromics-Bioinformatics/2.2_taxonomy/30_merge_summary/viral_contigs_all_outputs_merged.csv"
rafah_filename = "/home/qe76qox/Viromics-Bioinformatics/2.1_host_prediction/10_rafah/rafah_1_Host_Predictions_gtdb.csv"


# %%
df_jaeger = pd.read_csv(jaeger_filename, header=0, sep='\t',
                       usecols=[0,2,16]).set_index('contig_id')

df_jaeger = df_jaeger.loc[(df_jaeger['prediction'] == 'phage')]

print(df_jaeger)

# %%
df_assembly = pd.read_csv(assembly_filename, header=0, sep='\t',
                         usecols=[0,1,2,3],
                         names=['contig_id','length','coverage','is_circular']).set_index('contig_id')
print(df_assembly)

# %%
df_checkv = pd.read_csv(checkv_filename, header=0, sep='\t',
                       usecols=[0,4,6,7,9,11],
                       names=['contig_id','checkv_total_genes','checkv_viral_genes','checkv_quality','completeness','contamination']).set_index('contig_id')

df_checkv.fillna(0, inplace=True)
# print(df_checkv)

df_checkv = df_checkv.loc[(df_checkv['completeness'] > 50) & (df_checkv['contamination'] < 5)]

print(df_checkv)

# %%
df_pharokka = pd.read_csv(pharokka_filename, header=0, sep='\t')

df_pharokka = df_pharokka.pivot(index='contig',columns='Description')
df_pharokka.columns = df_pharokka.columns.get_level_values(1)

#df_pharokka = df_pharokka[['Count']]
df_pharokka = df_pharokka[['CDS','CRISPRs','DNA, RNA and nucleotide metabolism','connector', 'head and packaging','integration and excision', 'lysis','moron, auxiliary metabolic gene and host takeover','other','tRNAs','tail','transcription regulation']]

df_pharokka.rename(columns={'CDS':'phanotate_total_genes'}, inplace=True)


print(df_pharokka)

# %%
# genomad genes

df_genomad_genes = pd.read_csv(genomad_filename, header=0, sep='\t',
                      usecols=[0,8],
                      names=['contig_id_cds','marker'])
#print(df_genomad_genes)

# split the column, join to df and add a column name
df_genomad_genes = df_genomad_genes.join(df_genomad_genes['contig_id_cds'].str.rsplit('_',n=1, expand=True).add_prefix('A'))

# rename the columns
df_genomad_genes.rename(columns={'A0':'contig_id','A1':'cds'}, inplace=True)
# print(df_genomad_genes)

# get total genes
df_genomad_genes_counts = df_genomad_genes.groupby(by='contig_id').count()
# print(df_genomad_genes_counts)

# get viral genes base on the marker column containing "VV"
df_genomad_genes['marker'].fillna('bla',inplace=True)
df_genomad_vv_counts = df_genomad_genes[df_genomad_genes['marker'].str.contains('VV')].groupby(by='contig_id').count()
print(df_genomad_vv_counts)

# merge the gene counts and marker gene counts and rename
df_final_genomad = pd.DataFrame.join(df_genomad_genes_counts[['cds']],df_genomad_vv_counts[['marker']],on='contig_id', how='left')
df_final_genomad.rename(columns={'cds':'genomad_total_genes','marker':'genomad_viral_genes'}, inplace=True)

print(df_final_genomad)

# %%
df_rafah = pd.read_csv(rafah_filename, header=0,
                      usecols=[0,1,2],
                      names=['contig_id','rafah_predicted_host','rafah_prediction_score']).set_index('contig_id')
print(df_rafah)


# %%
# vcontact3 taxonomy

df_vcontact = pd.read_csv(vcontact_filename,header=0)
df_vcontact = df_vcontact[df_vcontact['GenomeName'].str.contains('contig')]

cols = ['realm (prediction)','phylum (prediction)','class (prediction)','order (prediction)','family (prediction)','subfamily (prediction)','genus (prediction)']
df_vcontact['taxonomy'] = df_vcontact[cols].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)

df_vcontact_final=df_vcontact[['GenomeName','taxonomy']]
df_vcontact_final.rename(columns={'GenomeName':'contig_id','taxonomy':'vcontact_taxonomy'}, inplace=True)
df_vcontact_final.set_index('contig_id',inplace=True)

print(df_vcontact_final)

# %%
# genomad taxonomy

df_genomad_taxa = pd.read_csv(genomad_taxonomy_filename,sep='\t',header=0,
                             usecols=[0,2,3,4], names=['contig_id','genomad_agreement','genomad_taxid','genomad_taxonomy']).set_index('contig_id')

print(df_genomad_taxa)

# %%
# merging all the files together

dfs_list = [df_assembly,df_pharokka,df_final_genomad, df_rafah,df_vcontact_final,df_genomad_taxa]
df_jaeger_checkv = df_checkv.join(df_jaeger, how='inner')
print(df_jaeger_checkv)


df_others = ft.reduce(lambda left, right: pd.merge(left,right,how='left',left_index=True,right_index=True), dfs_list)

df_final = df_jaeger_checkv.join(df_others,how='left')

print(df_final)

# %%
df_final.to_csv(out_filename)

