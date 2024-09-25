#!/bin/bash
#SBATCH --tasks=1
#SBATCH --job-name phold
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output 3_Project/3.1_phold/phold_%j_1.out
#SBATCH --error 3_Project/3.1_phold/phold_%j_1.err
#SBATCH --gres=gpu:1

# define directories
indir=1.4_annotation/10_pharokka
outdir=3_Project/3.1_phold

## load gpu module
module load nvidia/cuda/12.1.0

# conda environment
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh

conda activate pholdENV_v0.2.0_gpu

# have phold predict outputs based on the pharokka predictions
phold predict -i $indir/pharokka.gbk -o $outdir/10_phold_predict_output
# compare pholds output to the pharokka predictions
phold compare -i $indir/pharokka.gbk -o $outdir/20_phold_compare_output  --predictions_dir $outdir/10_phold_predict_output -t 32