#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=1G
#SBATCH --time=2:00:00
#SBATCH --job-name=sequencing_QC
#SBATCH --output=1.1_QC/10_nanoplot/nanoplot.slurm.out.%j
#SBATCH --error=1.1_QC/10_nanoplot/nanoplot.slurm.err.%j

# First, we check the quality of the reads using nanoplot
# activate the conda environment containing nanoplot

source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh && conda activate nanoplot_v1.41.3

datadir="./data/sequences"
outdir="./1.1_QC/10_nanoplot"

#check after filtering
#datadir="./1.1_QC/20_chopper"
#outdir="./1.1_QC/20_chopper/checked_files"

# Run nanoplot on 3 fastq files in a for loop
# -t : threads (cpus)
# --plots : types of plots to produce, see gallery: https://gigabaseorgigabyte.wordpress.com/2017/06/01/example-gallery-of-nanoplot/ 
# --N50 : draw the N50 vertical line on the plots
# --color : color for the plots
# -o : output directory
# to see other options for colors and prefilters, run : NanoPlot -h

for fn in $datadir/full_barcode6*.fastq.gz
do
base_f=$(basename $fn .fastq.gz)        # select just the first part of the name
NanoPlot -t 10 --plots dot --N50 --color slateblue --fastq $fn -o $outdir/${base_f}/
done