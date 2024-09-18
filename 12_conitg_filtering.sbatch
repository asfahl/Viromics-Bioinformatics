#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=short
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --job-name=filter_contigs
#SBATCH --output=1.3_virus_identification/40_results_filter_contigs/filter_contigs.slurm.%j.out
#SBATCH --error=1.3_virus_identification/40_results_filter_contigs/filter_contigs.slurm.%j.err

# activate the python virtual environment with the packages we need
source py3env/bin/activate

# in this sbatch script, its not necessarry to create the directory,
# we already told sbatch to create it for the log files.
# mkdir -p 40_results_filter_contigs

# In this solution, our script takes the assembly and the files
# 'complete_contigs_default_jaeger.tsv' from jaeger and 
# 'quality_summary.tsv' from CheckV as an input. Set them as variables
# for readability
cross_assembly='1.2_assembly/10_results_assembly_flye/cross_assembly/assembly.fasta'
jaegerresults='1.3_virus_identification/10_jaeger/results_jaeger/assembly/assembly_default_jaeger.tsv'
checkvresults='1.3_virus_identification/20_results_assessment_checkv/cross_assembly/quality_summary.tsv'
outdir='1.3_virus_identification/40_results_filter_contigs'

# run our script for filtering contigs based on the output of jaeger
# and CheckV as well as the output path. 
python python_scripts/virus_identification/contig_filtering.py $cross_assembly $jaegerresults $checkvresults $outdir/assembly.fasta

# deactivate the environment
deactivate