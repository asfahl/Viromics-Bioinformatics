#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=short,standard
#SBATCH --mem=5G
#SBATCH --time=01:00:00
#SBATCH --job-name=blastn
#SBATCH --output=./2.1_host_prediction/20_blastn/blastn.slurm.%j.out
#SBATCH --error=./2.1_host_prediction/20_blastn/blastn.slurm.%j.err

# reference to tool locations
makebalstdb=/home/groups/VEO/tools/blast/v2.14.0/bin/makeblastdb
blastn=/home/groups/VEO/tools/blast/v2.14.0/bin/blastn

# input files (metagenome bins and cross-assembly)
bins_fasta=2.1_host_prediction/20_blastn/good_bins_all.fa
assembly=1.3_virus_identification/40_results_filter_contigs/assembly.fasta

$makebalstdb -in $bins_fasta -out $bins_fasta -dbtype nucl

$blastn -db $bins_fasta -query $assembly -out viral_contigs_vs_bins_blastn.out -num_threads 32 -outfmt 6



