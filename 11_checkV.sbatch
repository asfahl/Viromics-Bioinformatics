#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=22
#SBATCH --partition=standard
#SBATCH --mem=20G
#SBATCH --time=02:30:00
#SBATCH --job-name=assessment_checkv
#SBATCH --output=1.3_virus_identification/20_results_assessment_checkv/assessment_checkv.slurm.%j.out
#SBATCH --error=1.3_virus_identification/20_results_assessment_checkv/assessment_checkv.slurm.%j.err



# activate the conda environment containing the checkv installation
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh && conda activate checkv_v1.0.1

# CheckV parameters (https://bitbucket.org/berkeleylab/checkv/src/master/#markdown-header-running-checkv)
# checkv end-to-end runs the CheckV pipeline from end to end :). It expects an input fasta file 
# with the assembly and an output path.
# -t: threads

# assigning variables for readability

database='/work/groups/VEO/databases/checkv/v1.5'
outdir='1.3_virus_identification/20_results_assessment_checkv'
cross_assembly='1.2_assembly/10_results_assembly_flye/cross_assembly/assembly.fasta'
#single_assembly='1.2_assembly/20_results_assessment_vclust/single_assemblies/assembly.fasta'

checkv end_to_end -t 20 -d $database $cross_assembly $outdir/cross_assembly
#checkv end_to_end -t 20 -d $database $single_assembly $outdir/single_assemblies
