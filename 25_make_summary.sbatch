#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=short
#SBATCH --mem=5G
#SBATCH --time=00:10:00
#SBATCH --job-name=merge
#SBATCH --output=2.2_taxonomy/30_merge_summary/merge.slurm.%j.out
#SBATCH --error=2.2_taxonomy/30_merge_summary/merge.slurm.%j.err

# activate the python virtual environment with the packages we need
source ./py3env/bin/activate


python3 python_scripts/taxonomy/merge_summary.py

# deactivate the environment
deactivate