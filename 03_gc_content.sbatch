#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --job-name=sequencing_QC
#SBATCH --output=1.1_QC/30_gc-content/gc-plot.slurm.out.%j
#SBATCH --error=1.1_QC/30_gc-content/gc-plot.slurm.err.%j

# activate the py3env
source ./py3env/bin/activate

# activate the py3env
source ./py3env/bin/activate

# first argument: ="./data/sequences"
# second argument: ="./data/sequences" 
python3 python_scripts/qc/gc_content.py ./data/sequences ./data/sequences/T4_SRR29341083.fastq.gz

deactivate