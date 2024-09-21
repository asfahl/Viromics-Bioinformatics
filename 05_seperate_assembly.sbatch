#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH --job-name=assembly_flye
#SBATCH --output=1.2_assembly/10_results_assembly_flye/single_assemblies/assembly_flye.slurm.%j.out
#SBATCH --error=1.2_assembly/10_results_assembly_flye/single_assemblies/assembly_flye.slurm.%j.err

# First, activate the conda environment which holds the flye installation on draco:
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh && conda activate flye_v2.9.2

# set a data directory holding your samples in .fastq.gz format
datadir="1.1_QC/20_chopper"

outdir="1.2_assembly/10_results_assembly_flye/single_assemblies"
mkdir -p $outdir
for barcode in $(seq 62 64)
do
        flye --nano-raw $datadir/full_barcode${barcode}_filtered_reads.fastq.gz --meta --out-dir $outdir/barcode$barcode -t 30
done

# close the conda environment
conda deactivate

# activate your virtual environment
source py3env/bin/activate

# rename the contings in the generated single assemblies with a python script.
# It assumes three parameters here, first the input file name, second the output file name and third,
# the sample name to add to each contig name in the respective fasta file.
for barcode in $(seq 62 64)
do
        python3 python_scripts/assembly/10_run_rename_scaffolds.py $outdir/barcode$barcode/assembly.fasta $outdir/barcode$barcode.fasta barcode$barcode
done

deactivate

cat $outdir/barcode*.fasta > $outdir/assembly.fasta