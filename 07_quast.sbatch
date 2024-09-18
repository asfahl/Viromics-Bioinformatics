#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --job-name=assessment_quast
#SBATCH --output=./1.2_assembly/30_results_assessment_quast/assessment_quast.slurm.%j.out
#SBATCH --error=./1.2_assembly/30_results_assessment_quast/assessment_quast.slurm.%j.err

# Set some variables for the quast script on draco and the files to be analysed
quast='python3 /home/groups/VEO/tools/quast/v5.2.0/quast.py'
cross_assembly='1.2_assembly/10_results_assembly_flye/cross_assembly/assembly.fasta'
single_assembly='1.2_assembly/20_results_assessment_vclust/single_assemblies/assembly.fasta'
outdir='1.2_assembly/30_results_assessment_quast'

# run Quast to visualize the distribution of contig lengths
$quast -o $outdir/cross_assembly $cross_assembly
$quast -o $outdir/single_assemblies $single_assembly

