#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=standard
#SBATCH --mem=5G
#SBATCH --time=02:00:00
#SBATCH --job-name=Jaeger
#SBATCH --output=1.3_virus_identification/10_jaeger/jaeger.slurm.out.%j
#SBATCH --error=1.3_virus_identification/10_jaeger/jaeger.slurm.err.%j

# set directories
indir="1.2_assembly/10_results_assembly_flye/cross_assembly/assembly.fasta"
outdir="1.3_virus_identification/10_jaeger"
mkdir -p $outdir


# activate cuda
module load nvidia/cuda/12.1.0

# activate jaeger
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh
conda activate jaeger_dev

# run jaeger
python3 /home/groups/VEO/tools/jaeger/v1.1.30a0/Jaeger/bin/jaeger run -i $indir -o $outdir/results_jaeger