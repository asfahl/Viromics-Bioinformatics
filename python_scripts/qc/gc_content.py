#! usr/bin/python3

# import libraries into your script
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
import pandas as pd
import numpy as np
import gzip, sys
import matplotlib.pyplot as plt

fastq_dir=sys.argv[1]
T4_path= sys.argv[2]

# read in fastq.gz files 
gc_62={}
with gzip.open(f"{fastq_dir}/full_barcode62.fastq.gz", 'rt') as f:
   for record in SeqIO.parse(f, "fastq"):
       gc_62[record.id]=gc_fraction(record.seq)*100

print("gc_62... done")

gc_63={}
with gzip.open(f"{fastq_dir}/full_barcode63.fastq.gz", 'rt') as f:
   for record in SeqIO.parse(f, "fastq"):
       gc_63[record.id]=gc_fraction(record.seq)*100
print("gc_63... done")

gc_64={}
with gzip.open(f"{fastq_dir}/full_barcode64.fastq.gz", 'rt') as f:
   for record in SeqIO.parse(f, "fastq"):
       gc_64[record.id]=gc_fraction(record.seq)*100
print("gc_64... done")

gc_T4={}
with gzip.open(T4_path, 'rt') as f:
   for record in SeqIO.parse(f, "fastq"):
       gc_T4[record.id]=gc_fraction(record.seq)*100
       #print(gc_fraction(record.seq)*100)
print("gc_T4... done")

# convert to a pandas dataframe
# For your samples
df_62 = pd.DataFrame([gc_62])
df_62 = df_62.T
df_63 = pd.DataFrame([gc_63])
df_63 = df_63.T
df_64 = pd.DataFrame([gc_64])
df_64 = df_64.T

# For T4 phage
df_T4 = pd.DataFrame([gc_T4])
df_T4 = df_T4.T

# make a plotting function
def make_plot(axs):
   # We can set the number of bins with the *bins* keyword argument.
   n_bins = 100

   ax1 = axs[0]
   ax1.hist(df_62, bins=n_bins)
   ax1.set_title('% GC content')
   ax1.set_ylabel('Barcode 62')
   ax1.set_xlim(0, 100)
   ax1.set_xticks(np.arange(0, 100, step=10))

   ax2 = axs[1]
   ax2.hist(df_63, bins=n_bins)
   ax2.set_ylabel('Barcode 63')

   ax3 = axs[2]
   ax3.hist(df_64, bins=n_bins)
   ax3.set_ylabel('Barcode 64')

   ax4 = axs[3]
   ax4.hist(df_T4, bins=n_bins)
   ax4.set_ylabel('T4 Phage')

# Plot:
fig, axs = plt.subplots(4,1, sharex=True, tight_layout=True)
make_plot(axs)
plt.show()

# save your plot
fig.savefig('GC_Content.png', dpi=150)





