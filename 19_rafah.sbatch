#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=standard
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name=rafah_host_prediction
#SBATCH --output=2.1_host_prediction/10_rafah/rafah.slurm.%j.out
#SBATCH --error=2.1_host_prediction/10_rafah/rafah.slurm.%j.err

assembly='1.3_virus_identification/40_results_filter_contigs/assembly.fasta'
contigs='2.1_host_prediction/10_rafah/split_contigs'

# activate py3env
source py3env/bin/activate

# run python script to seperate the assembly into individual genomes
python python_scripts/host_prediction/split_assembly_files.py $assembly $contigs

# deactivate py3env
deactivate

# activate the conda environment with the dependencies RaFAH requires
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh && conda activate perl_v5.32.1

# set a variable to call the RaFAH script
rafah='/home/groups/VEO/tools/rafah/RaFAH.pl'

# rafah parameters (https://gensoft.pasteur.fr/docs/RaFAH/0.3/)
# --predict: run the pipeline for predicting hosts
# --genomes_dir: the directory with the separate files for each contig
# --extension: the extension of the contig files
# --file_prefix: can specify an output dir here and "run name" (rafah_1) here
perl "$rafah" --predict --genomes_dir $contigs --extension .fasta --file_prefix 2.1_host_prediction/10_rafah/rafah_1

# deactivate conda environment
conda deactivate

# activate py3env again
source py3env/bin/activate

# translate the annotations from ncbi to gtdb
# set the path to the table with the translation between NCBI and GTDB and to RaFAH output
translation='2.1_host_prediction/10_rafah/ncbi_to_gtdb.csv'
rafahtable='2.1_host_prediction/10_rafah/rafah_1_Host_Predictions.tsv'

# run our script with the appropriate file names as parameters
python python_scripts/host_prediction/translate_taxonomy.py $rafahtable $translation 2.1_host_prediction/10_rafah/rafah_1_Host_Predictions_gtdb.csv

# deactivate py3env
deactivate