#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=standard,interactive
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --job-name=pharokka_genome_annotation
#SBATCH --output=3_Project/3.4_Dataset_30/10_pharokka/pharokka.slurm.%j.out
#SBATCH --error=3_Project/3.4_Dataset_30/10_pharokka/pharokka.slurm.%j.err

# set up directories
indir=3_Project/3.0_select_contigs
outdir=3_Project/3.4_Dataset_30/10_pharokka

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
pharokka.py -i $indir/select_contigs_30.fasta -o $outdir -d /work/groups/VEO/databases/pharokka/v1.7.1/ -t 10 -m --meta_hmm -f -g phanotate