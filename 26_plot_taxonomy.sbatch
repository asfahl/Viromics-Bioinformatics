#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=short
#SBATCH --mem=5G
#SBATCH --time=00:10:00
#SBATCH --job-name=taxonomic_graph
#SBATCH --output=2.2_taxonomy/40_plot_networks/taxonomic_graph.slurm.%j.out
#SBATCH --error=2.2_taxonomy/40_plot_networks/taxonomic_graph.slurm.%j.err

# activate the python virtual environment with the packages we need
source py3env/bin/activate

# in this sbatch script, its not necessarry to create the directory,
# we already told sbatch to create it for the log files.
mkdir -p 2.2_taxonomy/40_plot_networks/

# define an array with contig ids to plot
declare -a arr=("contig_871" "contig_267", "contig_825")

# loop through the array
for contig in "${arr[@]}"
do
	# run our scripts for plotting the network around a contig of choice
  # PLOT CONTIG IN LOCAL GENE SHARING NETWORK
	python python_scripts/visualization/taxonomy_visualization/plot_local_taxonomy.py 2.2_taxonomy/10_vcontact_results/graph.cyjs 2.2_taxonomy/40_plot_networks/"$contig"_graph_local.png $contig

  # PLOT CONTIG IN GLOBAL GENE SHARING NETWORK
  python python_scripts/visualization/taxonomy_visualization/plot_global_taxonomy.py 2.2_taxonomy/10_vcontact_results/graph.cyjs 2.2_taxonomy/40_plot_networks/"$contig"_graph_global.png $contig
done


# deactivate the environment
deactivate