#!/bin/bash
#SBATCH --tasks=1
#SBATCH --job-name phold
#SBATCH --partition=standard,interactive
#SBATCH --cpus-per-task=32
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output 3_Project/3.2_subsampled_phold/20_phold_plot/phold_%j_1.out
#SBATCH --error 3_Project/3.2_subsampled_phold/20_phold_plot/phold_%j_1.err

# define directories
indir=3_Project/3.2_subsampled_phold/10_phold_on_pharokka
outdir=3_Project/3.2_subsampled_phold/20_phold_plot

## load gpu module
module load nvidia/cuda/12.1.0

# conda environment
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh

conda activate pholdENV_v0.2.0_gpu

# replace phold predict and compare with phold run to do both for better optimisation on cpu
phold plot -i $indir/phold.gbk -o $outdir -f