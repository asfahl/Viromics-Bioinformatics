#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short
#SBATCH --mem=5G
#SBATCH --time=01:00:00
#SBATCH --job-name=combine
#SBATCH --output=./2.1_host_prediction/30_combined_results/combine.slurm.%j.out
#SBATCH --error=./2.1_host_prediction/30_combined_results/combine.slurm.%j.err

# Our python script for combining blastn and rafah needs our python environment again
source py3env/bin/activate

# file paths
blast='2.1_host_prediction/20_blastn/viral_contigs_vs_bins_blastn.out'
gtdb='2.1_host_prediction/20_blastn/gtdbtk.bac120.summary.tsv'
rafah='2.1_host_prediction/10_rafah/rafah_1_Host_Predictions_gtdb.csv'
output='2.1_host_prediction/30_combined_results/combined_blast_rafah_host_predictions.csv'

python3 python_scripts/host_prediction/merged_predictions.py $blast $gtdb $rafah $output
