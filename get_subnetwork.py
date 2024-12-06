import numpy as np
import pandas as pd
import networkx as nx
import sys
import os
from loopy_bp import get_prior_network

# functions
def spath(G, pair_list, max_path_length):
    #### get all shortest paths for a given gene list ####
    paths = []
    for pair in pair_list:
        # check for node, path and pathlength requirements:
        if G.has_node(pair[0]) and G.has_node(pair[1]):
            if nx.has_path(G, pair[0], pair[1]):
                if nx.shortest_path_length(G, pair[0], pair[1]) <= max_path_length:
                    # get all shortest paths for current pair
                    for path in nx.all_shortest_paths(G, pair[0], pair[1]):
                        # add each path as individual pairs so [0,1,2] becomes [0,1] and [1,2]
                        for i in range(len(path)-1):
                            paths.append([path[i], path[i+1]])
    
    # make each pair unique regardless of element order
    unique_paths = []
    for path in paths:
        same_path = 0
        for added_path in unique_paths:
            if set(path) == set(added_path):
                same_path = 1
        if same_path == 0:
            unique_paths.append(path)
    # output dataframe of edges for easy merge operations
    path_df = pd.DataFrame(unique_paths, columns=['Gene_A', 'Gene_B'])
    return path_df

if __name__=="__main__":
    # inputs:
    ppi_file = sys.argv[1]
    edge_file = sys.argv[2]
    out_path = sys.argv[3]
    prob_thresh = float(sys.argv[4])
    gene = sys.argv[5]
    max_path = 2
    
    # Step 1: subset edge data based on probability threshold
    # read in edge_data
    edge_data = pd.read_csv(edge_file)
    # subset edge data
    edge_data = edge_data.loc[edge_data["posterior_edge"]>=prob_thresh,]
    
    # Step 2: create G with remaing edges
    G = nx.from_pandas_edgelist(edge_data, source="Gene_A", target="Gene_B")
    
    # Step 3: use spath alg to get 2 hop expanded network from gene of interest
    # get unique non-identical pairs of gene list
    pair_list = [(gene, node) for node in G.nodes() if node != gene]
    path_AB = spath(G, pair_list, max_path)
    path_BA = pd.DataFrame({"Gene_A": path_AB["Gene_B"].tolist(), "Gene_B": path_AB["Gene_A"].tolist()})
    path_df = pd.concat([path_AB, path_BA])
    
    # Step 4: add back edge probability information
    final_edges = pd.merge(path_df, edge_data, how="left", on=["Gene_A", "Gene_B"])
    # remove non-existent edges using the NaNs to subset
    final_edges = final_edges.dropna()
    # gene specific output dir
    os.makedirs(f"{out_path}/20241030_Resutls/p-sites_of_interest", exist_ok=True)

    # save
    final_edges.to_csv(f"{out_path}/20241030_Resutls/p-sites_of_interest/{gene}_subnetwork.csv", index=False)
    