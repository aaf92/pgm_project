# function to make a gene cluster dictionary
# dictionary keys = cluster number, values = list of genes in cluster
make_clus_dict <- function(clus_path) {
  # initialize dictionary
  gene_clus_dict  <- list()
  
  # fill in with genes per cluster
  for (i in Sys.glob(clus_path)) {
    genes <- read.csv(i, header = FALSE, sep = '\t')
    clus_num <- as.integer(unlist(strsplit(unlist(strsplit(unlist(strsplit(i,"/"))[length(unlist(strsplit(i,"/")))],"clus"))[2],"_"))[1])
    gene_clus_dict[[clus_num]] <- genes[[1]]
  }
  return(gene_clus_dict)
}



# function to create an edge list from from site-to-gene weights
# sites = rows in weights variable, genes = cols in weights variables
make_edgelist <- function(weights) {
  # create edgelist
  edgelist <- which(weights != 0, arr.ind = TRUE)
  edgelist[,2] <- dim(weights)[1] + edgelist[,2]
  return(edgelist)
}



# function to get node colors for network
get_node_colors <- function(weights,gene_clus_dict,gene_colors,site_colors) {
  # get gene node colors
  clus_assign <- c()
  for (i in 1:dim(weights)[2]) {
    for (j in 1:length(gene_clus_dict)) {
      if (sum(gene_clus_dict[[j]] == colnames(weights)[i]) > 0) {
        clus_assign <- c(clus_assign,j)
        break
      } 
    }
  }
  node_colors <- c(rep(c(site_colors),each=dim(weights)[1]),clus_colors[clus_assign])
  return(node_colors)
}


