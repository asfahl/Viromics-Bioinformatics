import os, sys, json
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import numpy as np

# the keys to access taxonomic ranks of the phages in the graph
ranks = ["order", "class", "phylum"]

# this function looks for the lowest available rank according to the ones defined in the ranks variable
def get_lowest_tax(node):
    for r in ranks: 
        if r in node.keys(): return node[r]
    return "none"

def main():

    graph_filename = os.path.abspath(sys.argv[1])
    assert graph_filename.endswith(".cyjs")
    
    out_filename = os.path.abspath(sys.argv[2])
    assert out_filename.endswith(".png")

    contig_id = sys.argv[3]

    with open(graph_filename) as json_file:
        graph_json = json.load(json_file)

        # set for compatibility with networkx
        graph_json["data"] = []

        for node in graph_json["elements"]["nodes"]:
            # assign "value" attribute for compatibility
            node["data"]["value"] = node["data"]["id"]
            # assign lowest availble taxonomic rank with our function
            node["data"]["lowest_tax"] = get_lowest_tax(node["data"])

        g = nx.cytoscape_graph(graph_json)
        
    # find your favourite contig and get its internal node id
    contig_id_node = -1
    for node in g.nodes(data="name"):
    	if node[1] == contig_id: 
            contig_id_node = node[0]
    
    # raise an error, if the contig could not be found
    if contig_id_node == -1: raise(AssertionError(f"Could not find {contig_id} in the graph"))

    # get all nodes connected to the selected contig and compute the corresponding subgraph
    node_set = {contig_id_node}.union({n for n in g.neighbors(contig_id_node)})
    subgraph =  g.subgraph(node_set)

    # get all taxa in the subgraph
    taxa = nx.get_node_attributes(subgraph, "lowest_tax")
    
    # set up a mapping from node to color based on the taxonomic rank
    unique_taxa = np.unique([t for t in taxa.values()])
    taxa_to_color = {t:mpl.colormaps["tab20"](i) for i, t in enumerate(unique_taxa)}
    colors = {n:taxa_to_color[t] for n, t in taxa.items()}

    # compute a positional layout using the kamada-kawai algorithm
    pos = nx.kamada_kawai_layout(subgraph)

    # set colors for nodes. use the taxonomic information to color the reference nodes
    names = nx.get_node_attributes(subgraph, "name")
    colors = [('seagreen' if names[n] == contig_id else 'red') if names[n].startswith('contig') else colors[n] for n in subgraph.nodes()]

    # draw the network using the functions networkx provides
    nx.draw_networkx_edges(subgraph, pos=pos, width=0.05, alpha=0.4)
    nx.draw_networkx_nodes(subgraph, pos=pos, alpha=0.7, node_size=100, node_color=colors)
    nx.draw_networkx_labels(subgraph, pos=pos, labels=names, font_size=3)

    # create rectangles for a legend connecting taxa with their colors in the graph
    patches = [mpatches.Patch(color=c, label=l, alpha=0.6) for l,c in taxa_to_color.items() if l!='none'] + \
   	   [mpatches.Patch(color='red', label='our contigs'), mpatches.Patch(color='green', label='selected contig')]
   	
    plt.legend(handles=patches, loc='upper right', framealpha=0.5, frameon=True, prop={'size': 5})

    ax = plt.gca()
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(out_filename, dpi=700)


if __name__ == "__main__":
    main()

