#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
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

source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
conda activate pholdENV
echo "activated env"

phold run -i $indir/contig_871.gbk --output_dir=$outdir -t 10