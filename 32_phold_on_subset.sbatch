#!/bin/bash
#SBATCH --tasks=1
#SBATCH --job-name phold
#SBATCH --partition=standard,interactive
#SBATCH --cpus-per-task=32
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output 3_Project/3.4_Dataset_30/20_phold_pharokka/phold_%j_1.out
#SBATCH --error 3_Project/3.4_Dataset_30/20_phold_pharokka/phold_%j_1.err

# define directories
indir=3_Project/3.4_Dataset_30/10_pharokka
outdir=3_Project/3.4_Dataset_30/20_phold_pharokka

## load gpu module
module load nvidia/cuda/12.1.0

# conda environment
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh

conda activate pholdENV_v0.2.0_gpu

# replace phold predict and compare with phold run to do both for better optimisation on cpu
phold run -d /veodata/01/databases/phold/v240724/ -i $indir/pharokka.gbk -o $outdir -t 96 --cpu -f