#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=short,standard,interactive
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --job-name=alignment_minimap2
#SBATCH --output=./1.2_assembly/40_results_alignment_minimap2/alignment_minimap2.slurm.%j.out
#SBATCH --error=./1.2_assembly/40_results_alignment_minimap2/alignment_minimap2.slurm.%j.err


# assign tool paths to aliases for better readability
minimap2='/home/groups/VEO/tools/minimap2/v2.26/minimap2'
samtools='/home/groups/VEO/tools/samtools/v1.17/bin/samtools'

# run minimap2 to align all reads from a sample to the assembled contigs
# and pipe the output into samtools for conversion into the binary bam format
#
# minimap2 parameters (https://lh3.github.io/minimap2/minimap2.html):
# -x map_ont : Use a preset for parameterizing the affine gap penalty model for the extension of matched seeds
# suited for noisy nanopore reads.
# -a : output in SAM format
# -t 30 : run with 30 threads
#
# samtools parameters (http://www.htslib.org/doc/samtools.html):
#
# samtools view can be used to convert between SAM, BAM and CRAM formats.
# view -u : output uncompressed binary format (BAM)
#
# samtools sort can be used to sort a SAM, BAM or CRAM file. Some tools expect sorted alignments.
# sort --write-index : output the index of the sorted alignments, can reduce file IO when accessing only a subset of the alignments
# sort -o : set the output file for the sorted alignments
#
# - : the - tells samtools to take the inpute from the pipe (| is the piping operator).


indir='1.1_QC/20_chopper'

# cross-assembly
outdir='data/alignments/cross_assembly'
assembly='1.2_assembly/10_results_assembly_flye/cross_assembly/assembly.fasta'
mkdir -p $outdir

# loop through the numbers 62 to 64 and use it to generate diffenrent filenames within the loop
for barcode in $(seq 62 64) 
do 
 $minimap2 -x map-ont -a -t 30  $assembly $indir/barcode${barcode}_filtered_reads.fastq.gz | \
   $samtools view -u - | $samtools sort -o $outdir/barcode$barcode.bam --write-index -
 $samtools stats $outdir/barcode$barcode.bam > $outdir/barcode$barcode_stats.txt
done

# single assemblies
outdir='data/alignments/single_assemblies'
assembly='1.2_assembly/20_results_assessment_vclust/single_assemblies/assembly.fasta'
mkdir -p $outdir

for barcode in $(seq 62 64) 
do 
 $minimap2 -x map-ont -a -t 30 $assembly $indir/barcode${barcode}_filtered_reads.fastq.gz | \
   $samtools view -u - | $samtools sort -o $outdir/barcode$barcode.bam --write-index -
 $samtools stats $outdir/barcode$barcode.bam > $outdir/barcode$barcode_stats.txt
done

