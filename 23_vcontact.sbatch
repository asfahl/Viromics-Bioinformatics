#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --job-name=vcontact3
#SBATCH --output=2.2_taxonomy/10_vcontact_results/vcontact3.slurm.out.%j
#SBATCH --error=2.2_taxonomy/10_vcontact_results/vcontact3.slurm.err.%j


# activate conda environment
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh && conda activate mamba_20231101_python_3.9

database='/work/groups/VEO/databases/vcontact3/v20231101/'
assembly='1.3_virus_identification/40_results_filter_contigs/assembly.fasta'
outdir='2.2_taxonomy/10_vcontact_results'

#options
# -e : how do you want the output formatted? See the vcontact3 help for more options

vcontact3 run --nucleotide $assembly --output $outdir --pyrodigal-gv --db-path $database -e cytoscape 