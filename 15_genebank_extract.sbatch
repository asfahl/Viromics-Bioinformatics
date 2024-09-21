#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --job-name=assembly_visualization
#SBATCH --output=1.4_annotation/10_pharokka/11_selected_contig/genebank_extraction.slurm.%j.out
#SBATCH --error=1.4_annotation/10_pharokka/11_selected_contig/genebank_extraction.slurm.%j.err

# set directories
indir=1.4_annotation/10_pharokka
outdir=1.4_annotation/10_pharokka/11_selected_contig

# activate the python virtual environment with the packages we need
source ./py3env/bin/activate

# change the output and contig names accordingly
python3 python_scripts/gene_annotation/gene_bank_extraction.py 1.4_annotation/10_pharokka/pharokka.gbk 1.4_annotation/10_pharokka/11_selected_contig/contig_871.gbk contig_871