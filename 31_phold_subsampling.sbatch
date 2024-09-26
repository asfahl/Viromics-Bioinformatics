#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=2G
#SBATCH --time=00:02:00
#SBATCH --job-name=py_script
#SBATCH --output=3_Project/3.0_select_contigs/extract_contigs.slurm.out.%j
#SBATCH --error=3_Project/3.0_select_contigs/extract_contigs.slurm.err.%j

# activate the python virtual environment with the packages we need
source py3env/bin/activate

# change the output and contig names accordingly
python3 python_scripts/Project/subsampling.py 1.4_annotation/10_pharokka/pharokka.gbk 3_Project/3.0_select_contigs/select_contigs_30.gbk