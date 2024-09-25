#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --job-name=genomad
#SBATCH --output=2.2_taxonomy/20_genomad_results/genomad.slurm.out.%j
#SBATCH --error=2.2_taxonomy/20_genomad_results/genomad.slurm.err.%j

# activate conda environment
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh && conda activate genNomad_v20230721

database='/work/groups/VEO/databases/geNomad/v1.3'
assembly='./1.3_virus_identification/40_results_filter_contigs/assembly.fasta'
outdir='./2.2_taxonomy/30_genomad_results'

#options

genomad annotate --threads 10 $assembly $outdir $database

