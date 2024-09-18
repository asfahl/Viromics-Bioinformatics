#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=1G
#SBATCH --time=2:00:00
#SBATCH --job-name=sequencing_QC
#SBATCH --output=1.1_QC/10_nanoplot/nanoplot.slurm.out.%j
#SBATCH --error=1.1_QC/10_nanoplot/nanoplot.slurm.err.%j

# activate your virtual environment
source py3env/bin/activate

python3 python_scripts/virus_identification/jaeger_analysis.py
deactivate