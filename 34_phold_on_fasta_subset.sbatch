#!/bin/bash
#SBATCH --tasks=1
#SBATCH --job-name phold
#SBATCH --partition=standard,interactive
#SBATCH --cpus-per-task=32
#SBATCH --mem=50G
#SBATCH --time=03:00:00
#SBATCH --output 3_Project/3.4_Dataset_30/30_phold_fasta/phold_fasta_%j_1.out
#SBATCH --error 3_Project/3.4_Dataset_30/30_phold_fastaphold_fasta_%j_1.err

# define directories
indir=3_Project/3.0_select_contigs
outdir=3_Project/3.4_Dataset_30/30_phold_fasta

## load gpu module
module load nvidia/cuda/12.1.0

# conda environment
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh

conda activate pholdENV_v0.2.0_gpu

# replace phold predict and compare with phold run to do both for better optimisation on cpu
phold run -d /veodata/01/databases/phold/v240724/ -i $indir/select_contigs_30.fasta -o $outdir -t 96 --cpu -f