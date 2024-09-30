#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=standard,interactive
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --job-name=phold_genome_annotation
#SBATCH --output=1.4_annotation/10_pharokka/12_selected_contig_phold/phold.slurm.%j.out
#SBATCH --error=1.4_annotation/10_pharokka/12_selected_contig_phold/phold.slurm.%j.err

# set up directories
indir=1.4_annotation/10_pharokka/11_selected_contig
outdir=1.4_annotation/10_pharokka/12_selected_contig_phold

# set up the environment
# conda info --envs: ich finde fold nicht?
deactivate

## load gpu module
module load nvidia/cuda/12.1.0

# conda environment
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh

conda activate pholdENV_v0.2.0_gpu

# replace phold predict and compare with phold run to do both for better optimisation on cpu
phold run -d /veodata/01/databases/phold/v240724/ -i $indir/pharokka.gbk -o $outdir -t 96 --cpu -f

# replace phold predict and compare with phold run to do both for better optimisation on cpu
#phold plot -i $outdir/phold.gbk -o $outdir/plot