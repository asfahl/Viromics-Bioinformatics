#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=standard,interactive
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --job-name=pharokka_genome_annotation
#SBATCH --output=1.4_annotation/10_pharokka/pharokka.slurm.%j.out
#SBATCH --error=1.4_annotation/10_pharokka/pharokka.slurm.%j.err

# set up directories
indir=1.3_virus_identification/40_results_filter_contigs
outdir=1.4_annotation/10_pharokka

# set up the environment
deactivate

source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
conda activate pharokka_v1.7.1

# run pharokka
# -i input
# -o output
# -d database directory
# -g use phanota for metagenom
# --meta_hmm sets the HMM
pharokka.py -i $indir/assembly.fasta -o $outdir -d /work/groups/VEO/databases/pharokka/v1.7.1/ -t 10 -m --meta_hmm -f -g phanotate