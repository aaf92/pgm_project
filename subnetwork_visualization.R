#### network code
setwd("/ix/djishnu/Aaron_F/PGM_project")

# libraries 
library("qgraph")
source("qgraph_funcs.R")

# paths
edge_path = "/ix/djishnu/Aaron_F/PGM_project/20241030_Resutls/p-sites_of_interest/Pou5f1_subnetwork.csv"

# load edge data
edge_data = read.csv(edge_path)


# network setup
# nodes
nodes <- unique(c(edge_data$Gene_A, edge_data$Gene_B))
node_colors <- rep("green", length(nodes))  # Set all nodes to blue
input_node <- "Pou5f1"        # Define the specific node
node_colors[which(nodes == input_node)] <- "red"  # Set specific node to red

# adjacency matrix for edges
adj_matrix <- matrix(0, nrow=length(nodes), ncol=length(nodes), dimnames=list(nodes, nodes))

# Populate the adjacency matrix based on directed edges
for (i in 1:nrow(edge_data)) {
  source_node <- as.character(edge_data$Gene_A[i])
  target_node <- as.character(edge_data$Gene_B[i])
  adj_matrix[source_node, target_node] <- 1  # Set 1 for a directed edge
}

edge_color <- '#dbd9d9'


# plot network
pdf("/ix/djishnu/Aaron_F/PGM_project/20241030_Resutls/p-sites_of_interest/Pou5f1_subnetwork_viz.pdf")
qgraph(adj_matrix, layout = "spring", directed = TRUE, color = node_colors, shape = "circle", edge.color = edge_color, label.scale = FALSE, label.cex = 0.5)
dev.off()

# pdf("20240709_pls1-3_avg_mod_neg_site_hit_weights_75perc_13modules.pdf")
# qgraph(neg_edgelist, layout = "spring", directed = FALSE, labels = neg_node_labels, color = neg_node_colors, 
#        shape = neg_node_shapes, node.width = neg_node_sizes, node.height = neg_node_sizes, edge.color = edge_color, 
#        label.scale = FALSE, label.cex = 0.5)
# dev.off()


