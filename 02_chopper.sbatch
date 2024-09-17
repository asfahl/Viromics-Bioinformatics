#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=1G
#SBATCH --time=2:00:00
#SBATCH --job-name=sequencing_QC
#SBATCH --output=1.1_QC/20_chopper/chopper.slurm.out.%j
#SBATCH --error=1.1_QC/20_chopper/chopper.slurm.err.%j

# activate chopper in the anaconda environment
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh && conda activate chopper_v0.5.0

# set directories
datadir="./data/sequences"
outdir="./1.1_QC/20_chopper"

# run chopper on a loop over the 3 fastq files
#   -q read quality
#   -l minimum length 
#   --maxlength
for fn in ./data/sequences/*.fastq.gz
do

base_f=$(basename $fn .fastq.gz)        # select just the first part of the name
gunzip -c $fn | chopper -q 14 -l 1000 --maxlength 40000 | gzip > $outdir/${base_f}_filtered_reads.fastq.gz
done

