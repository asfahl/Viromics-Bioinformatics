#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --job-name=assembly_visualization
#SBATCH --output=python_scripts/visualization/virome_metagenome.slurm.%j.out
#SBATCH --error=python_scripts/visualization/virome_metagenome.slurm.%j.err

indir=1.3_virus_identification/40_results_filter_contigs
refdir=python_scripts/visualization

source py3env/bin/activate

python3 python_scripts/visualization/virus_id.py $indir/assembly.fasta $refdir/input/metagenome.fasta $refdir/output/assembly_difference.png

deactivate

