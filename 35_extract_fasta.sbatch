#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=2G
#SBATCH --time=00:02:00
#SBATCH --job-name=py_script
#SBATCH --output=3_Project/3.0_select_contigs/extract_contigs_fasta.slurm.out.%j
#SBATCH --error=3_Project/3.0_select_contigs/extract_contigs_fasta.slurm.err.%j

declare -a arr=("contig_114" "contig_1393" "contig_1160" "contig_103" "contig_185" "contig_1169" "contig_200" "contig_92" "contig_1460" "contig_1552" "contig_1439" "contig_1323" "contig_1433" "contig_1033" "contig_1554" "contig_1547" "contig_58" "contig_256" "contig_601" "contig_825" "contig_1096" "contig_644" "contig_1131" "contig_1298" "contig_669" "contig_1545" "contig_815" "contig_1062" "contig_1170" "contig_267")

for conitg in "${arr[@]}"
do
  sed -n '/\<'$conitg'\>/,/>/p' 1.3_virus_identification/40_results_filter_contigs/assembly.fasta | head -n -1 >> 3_Project/3.0_select_contigs/select_contigs_30.fasta
  echo $conitg
done
