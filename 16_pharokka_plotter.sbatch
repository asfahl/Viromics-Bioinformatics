#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --job-name=assembly_visualization
#SBATCH --output=1.4_annotation/10_pharokka/11_selected_contig/pharokka_ploter.slurm.%j.out
#SBATCH --error=1.4_annotation/10_pharokka/11_selected_contig/pharokka_ploter.slurm.%j.err

# directories
indir=1.4_annotation/10_pharokka/11_selected_contig
outdir=1.4_annotation/10_pharokka/11_selected_contig

# activate the environment
# set up the environment
deactivate

source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
conda activate pharokka_v1.7.1

pharokka_plotter.py -i $indir/contig_871.fasta -n $outdir/contig_871_plot --gff $indir/contig_871.gff --genbank $indir/contig_871.gbk


