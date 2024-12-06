library(KEGGgraph)


kgml <- parseKGML("/Users/aaronfrancis/Downloads/mmu04620.xml")
nodes <- nodes(kgml)
edges <- edges(kgml)

edge_type <- c()
edge_subtype <- c()
edge_val <- c()
source <- c()
target <- c()
for (i in seq_along(edges)) {
  cur_type <- getType(edges[[i]])
  cur_subtype <- getSubtype(edges[[i]])
  
  # store edge_type
  edge_type <- c(edge_type, cur_type)
  
  # do a check for cur_subtype and store subtype name and value info
  if (is.null(cur_subtype) | length(cur_subtype) == 0) {
    edge_subtype <- c(edge_subtype, "")
    edge_val <- c(edge_val, "")
  }
  else {
    edge_subtype <- c(edge_subtype, getName(cur_subtype[[1]]))
    edge_val <- c(edge_val, getValue(cur_subtype[[1]]))
  }
  
  # store source and target node entry ids
  source <- c(source, getEntryID(edges[[i]])[[1]])
  target <- c(target, getEntryID(edges[[i]])[[2]])
}

source_names <- c()
target_names <- c()
for (i in seq_along(target)) {
  source_name <- strsplit(getDisplayName(nodes[source[i]][[1]]), ",")[[1]][1]
  target_name <- strsplit(getDisplayName(nodes[target[i]][[1]]), ",")[[1]][1]
  source_names <- c(source_names, source_name)
  target_names <- c(target_names, target_name)
}

# Create a dataframe without column names
pathway_net <- data.frame(source_names, target_names, edge_type, edge_subtype, edge_val)

# Specify column names
colnames(pathway_net) <- c("source", "target", "edge_type", "edge_name", "edge_value")

# drop nas
pathway_net_final <- na.omit(pathway_net)

# save
write.csv(pathway_net_final, "tlr_signaling_parsed_network_kegg_mmu04620.csv", row.names = FALSE)
