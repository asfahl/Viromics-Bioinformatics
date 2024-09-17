#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH --job-name=assembly_flye
#SBATCH --output=1.2_assembly/10_results_assembly_flye/cross_assembly/assembly_flye.slurm.%j.out
#SBATCH --error=1.2_assembly/10_results_assembly_flye/cross_assembly/assembly_flye.slurm.%j.err



# activate conda environment with flye installation on Draco
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh && conda activate flye_v2.9.2

# set directories
datadir="/home/qe76qox/Viromics-Bioinformatics/1.1_QC/20_chopper"
outdir="/home/qe76qox/Viromics-Bioinformatics/1.2_assembly/10_results_assembly_flye/cross_assembly"

# merge the sequences
cat $datadir/*.fastq.gz > $datadir/all_samples.fastq.gz

# create a folder for the output cross-assembly
mkdir -p $outdir


# complete the flye command. This is the computationally expensive part
# and profits from many cores (30 is a good number). The used memory should
# not exceed 20GB of RAM.
flye --nano-raw $datadir/all_samples.fastq.gz --meta --out-dir $outdir -t 30