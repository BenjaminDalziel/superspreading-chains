library(igraph)

# Taube, J.C., Miller, P.B. and Drake, J.M. (2022).
# An open-access database of infectious disease transmission trees
# to explore superspreader epidemiology. PLoS Biology, 20(6), e3001685.
# https://doi.org/10.1371/journal.pbio.3001685
# https://outbreaktrees.ecology.uga.edu/
# https://github.com/DrakeLab/taube-transmission-trees

# Read in data
tree_data <- readRDS("data/taube_2020_transmission_trees.rds")



# Extract metadata
metadata <- data.frame(
  name = tree_data$id,
  size = tree_data$Size,
  generations = tree_data$Generations,
  disease = tree_data$Disease,
  year = tree_data$Year,
  country = tree_data$Country,
  taube_R0 = tree_data$R0,
  taube_dispersion = tree_data$DP,
  source = tree_data$Source
)
ntree <- nrow(metadata)



# Storage
trees <- vector("list", ntree)



# Load tree(s)
for (i in 1:ntree) {
  # Load the tree
  # which is in a old igraph format
  # and gives a warning message
  # each time it is called.
  #
  # We suppress the messages
  # and take the tree to its adjacency matrix
  #
  # Then rebuild it such that
  # vertex names are their node numbers

  suppressMessages({
    this_tree <- tree_data[which(tree_data$id == metadata$name[i]), ]$tree[[1]]
    adj_mat <- as_adjacency_matrix(this_tree, sparse = FALSE)
  })

  row.names(adj_mat) <- seq_len(nrow(adj_mat))
  colnames(adj_mat) <- seq_len(ncol(adj_mat))
  this_tree <- graph_from_adjacency_matrix(adj_mat, mode = "directed")


  # Find out degrees and parent out degrees
  n <- nrow(adj_mat)
  out_degree <- .rowSums(adj_mat, n, n)
  z <- out_degree
  zp <- rep(NA, n)
  p <- rep(NA, n)
  for (j in 1:n) {
    pj <- which(adj_mat[, j] == 1)
    if (length(pj)) {
      p[j] <- pj
      zp[j] <- z[p[j]]
    }
  }
  out_degree_parent <- zp
  index_parent <- p


  # Store results
  trees[[i]] <- list(
    tree_name = metadata$name[i],
    tree = this_tree,
    out_degree = out_degree,
    out_degree_parent = out_degree_parent,
    index_parent = p
  )
}

# Save
file <- "data/trees_processed.Rdata"
save(trees, metadata, file = file)
