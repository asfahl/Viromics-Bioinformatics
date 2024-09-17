#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --job-name=assessment_vclust
#SBATCH --output=./1.2_assembly/20_results_assessment_vclust/assessment_vclust.slurm.%j.out
#SBATCH --error=./1.2_assembly/20_results_assessment_vclust/assessment_vclust.slurm.%j.err


vclust='python3.9 /home/groups/VEO/tools/vclust/v1.0.3/vclust.py'

indir='1.2_assembly/10_results_assembly_flye/cross_assembly'
outdir='1.2_assembly/20_results_assessment_vclust/cross_assembly'
mkdir -p $outdir

$vclust prefilter -i $indir/assembly.fasta -o $outdir/fltr.txt
$vclust align -i $indir/assembly.fasta -o $outdir/ani.tsv -t 30 --filter $outdir/fltr.txt
$vclust cluster -i $outdir/ani.tsv -o $outdir/clusters_tani_90.tsv --ids $outdir/ani.ids.tsv --metric tani --tani 0.90
$vclust cluster -i $outdir/ani.tsv -o $outdir/clusters_tani_70.tsv --ids $outdir/ani.ids.tsv --metric tani --tani 0.70
$vclust cluster -i $outdir/ani.tsv -o $outdir/clusters_ani_90.tsv --ids $outdir/ani.ids.tsv --metric ani --ani 0.90
$vclust cluster -i $outdir/ani.tsv -o $outdir/clusterreps.tsv --ids $outdir/ani.ids.tsv --algorithm uclust --metric ani --ani 0.95 --cov 0.85 --out-repr

source py3env/bin/activate

python3 python_scripts/assembly/20_run_filter_representatives.py $indir/assembly.fasta $outdir/assembly.fasta $outdir/clusterreps.tsv

deactivate

echo Done running python script for cross-assembly...

indir='1.2_assembly/10_results_assembly_flye/single_assemblies'
outdir='1.2_assembly/20_results_assessment_vclust/single_assemblies'
mkdir -p $outdir

$vclust prefilter -i $indir/assembly.fasta -o $outdir/fltr.txt
$vclust align -i $indir/assembly.fasta -o $outdir/ani.tsv -t 30 --filter $outdir/fltr.txt
$vclust cluster -i $outdir/ani.tsv -o $outdir/clusters_tani_90.tsv --ids $outdir/ani.ids.tsv --metric tani --tani 0.90
$vclust cluster -i $outdir/ani.tsv -o $outdir/clusters_tani_70.tsv --ids $outdir/ani.ids.tsv --metric tani --tani 0.70
$vclust cluster -i $outdir/ani.tsv -o $outdir/clusters_ani_90.tsv --ids $outdir/ani.ids.tsv --metric ani --ani 0.90
$vclust cluster -i $outdir/ani.tsv -o $outdir/clusterreps.tsv --ids $outdir/ani.ids.tsv --algorithm uclust --metric ani --ani 0.95 --cov 0.85 --out-repr

source py3env/bin/activate

python3 python_scripts/assembly/20_run_filter_representatives.py $indir/assembly.fasta $outdir/assembly.fasta $outdir/clusterreps.tsv

deactivate

