library(igraph)

# Parameters
seed <- 123
nrandom <- 10000
min_generations <- 2
min_size <- 20
cap <- 21
parent_cutoffs <- 0:(cap - 1)
child_cutoffs <- c(0, 1, 4)
y_max <- c(0.3, 0.2, 0.1)




# Create a rooted tree from a Prufer sequence
make_rooted_from_prufer <- function(p, root, shuffle = FALSE) {
  if (shuffle) p <- sample(p)
  g <- make_from_prufer(p)
  parent <- bfs(g, root = root, parent = TRUE)$parent
  child <- setdiff(seq_len(vcount(g)), root)
  graph_from_edgelist(cbind(parent[child], child), directed = TRUE)
}


# Count dyads by parent (row) and child (column) out-degree, indexed from zero
dyad_counts <- function(g) {
  z <- degree(g, mode = "out")
  e <- ends(g, E(g), names = FALSE)

  width <- max(c(0L, z)) + 1
  index <- z[e[, 1]] + 1 + width * z[e[, 2]]
  matrix(tabulate(index, nbins = width^2), width, width)
}


# Bin degrees at or above cap into the final row and column
bin_degrees <- function(dyads, cap) {
  labels <- c(0:(cap - 1), paste0(cap, "+"))
  binned <- matrix(0, cap + 1, cap + 1,
    dimnames = list(parent = labels, child = labels)
  )
  for (n in 0:(nrow(dyads) - 1)) {
    for (m in 0:(ncol(dyads) - 1)) {
      i <- min(n, cap) + 1
      j <- min(m, cap) + 1
      binned[i, j] <- binned[i, j] + dyads[n + 1, m + 1]
    }
  }
  binned
}


# Count children whose own and parent degrees exceed the given cutoffs
count_children <- function(dyads, child_cutoff, parent_cutoffs) {
  parent_degree <- seq_len(nrow(dyads)) - 1
  child_degree <- seq_len(ncol(dyads)) - 1
  dyads <- dyads[, child_degree > child_cutoff, drop = FALSE]
  vapply(parent_cutoffs, function(cutoff) {
    sum(dyads[parent_degree > cutoff, , drop = FALSE])
  }, numeric(1))
}


# Load trees meeting the generation and size criteria
load("data/trees_processed.Rdata")
keep <- metadata$generations >= min_generations & metadata$size >= min_size
trees <- trees[keep]


# Prepare each tree's dyad counts, root, and Prufer sequence
tree_dyads <- prufer <- vector("list", length(trees))
roots <- integer(length(trees))
for (i in seq_along(trees)) {
  g <- upgrade_graph(trees[[i]]$tree)
  tree_dyads[[i]] <- dyad_counts(g)
  roots[i] <- which(degree(g, mode = "in") == 0)
  prufer[[i]] <- to_prufer(as_undirected(g))
}


# Pool observed dyad counts
width <- max(vapply(tree_dyads, nrow, integer(1)))
observed_dyads <- matrix(0, width, width)
for (dyads in tree_dyads) {
  idx <- seq_len(nrow(dyads))
  observed_dyads[idx, idx] <- observed_dyads[idx, idx] + dyads
}


# Count qualifying children in the observed trees
observed_counts <- vapply(child_cutoffs, count_children, numeric(cap),
  dyads = observed_dyads, parent_cutoffs = parent_cutoffs
)
# Include every child degree
total_counts <- count_children(observed_dyads, -1, parent_cutoffs)


# Count qualifying children in randomized trees
null_counts <- array(0, dim = c(cap, nrandom, length(child_cutoffs)))
set.seed(seed)
for (k in seq_len(nrandom)) {
  cat("\rNull replicate ", k, "/", nrandom, sep = "")
  flush.console()
  null_dyads <- matrix(0, width, width)

  for (i in seq_along(trees)) {
    null_tree <- make_rooted_from_prufer(prufer[[i]], roots[i], shuffle = TRUE)
    dyads <- dyad_counts(null_tree)
    idx <- seq_len(nrow(dyads))
    null_dyads[idx, idx] <- null_dyads[idx, idx] + dyads
  }

  null_counts[, k, ] <- vapply(child_cutoffs, count_children, numeric(cap),
    dyads = null_dyads, parent_cutoffs = parent_cutoffs
  )
}
cat("\n")


# Calculate observed and null proportions among qualifying parents
observed_prob <- observed_counts / total_counts
null_prob <- null_counts / total_counts
null_mean <- apply(null_prob, c(1, 3), mean)
null_interval <- apply(
  null_prob, c(1, 3), quantile, c(0.025, 0.1, 0.9, 0.975)
)


# Prepare binned dyad counts and heatmap labels
binned_dyads <- bin_degrees(observed_dyads, cap)
cells <- which(binned_dyads > 0, arr.ind = TRUE)
counts <- binned_dyads[cells]
label_colors <- ifelse(
  log1p(counts) > 0.6 * max(log1p(binned_dyads)),
  "white", "black"
)


# Tabulate trees, parents, children at each parent-degree threshold
parent_degree <- seq_len(width - 1)
parents_by_degree <- rowSums(observed_dyads)[-1] / parent_degree
tree_max_degree <- vapply(tree_dyads, nrow, integer(1)) - 1
count_table <- data.frame(
  parent_degree_threshold = parent_cutoffs, qualifying_parents = 0,
  qualifying_trees = 0, children_of_qualifying_parents = total_counts
)
for (i in seq_along(parent_cutoffs)) {
  count_table$qualifying_parents[i] <- sum(
    parents_by_degree[parent_degree > parent_cutoffs[i]]
  )
  count_table$qualifying_trees[i] <- sum(tree_max_degree > parent_cutoffs[i])
}
for (j in seq_along(child_cutoffs)) {
  count_table[[paste0("children_with_out_degree_gt_", child_cutoffs[j])]] <-
    observed_counts[, j]
}

dir.create("tables", showWarnings = FALSE)
write.csv(count_table, "tables/count_table.csv", row.names = FALSE)


# Plot observed proportions against null means and intervals
par(mar = c(3.8, 4.6, 1, 0.6),
  cex = 1, cex.axis = 1.1, cex.lab = 1.2, mgp = c(2, 0.7, 0)
)
for (j in seq_along(child_cutoffs)) {
  plot(
    parent_cutoffs, observed_prob[, j],
    type = "n",
    xlim = c(0, 20),
    ylim = c(0, y_max[j]),
    bty = "l",
    xaxt = "n",
    yaxt = "n",
    xlab = expression(z[p] > x),
    ylab = ""
  )
  axis(1, at = c(0, 5, 10, 15, 20))
  axis(2, at = c(0, y_max[j] / 2, y_max[j]), las = 1)
  title(ylab = bquote(P(z > .(child_cutoffs[j]))), line = 2.5)
  polygon(
    c(parent_cutoffs, rev(parent_cutoffs)),
    c(null_interval[1, , j], rev(null_interval[4, , j])),
    col = "#D0D0D0", border = NA
  )
  polygon(
    c(parent_cutoffs, rev(parent_cutoffs)),
    c(null_interval[2, , j], rev(null_interval[3, , j])),
    col = "#777777", border = NA
  )
  #lines(parent_cutoffs, null_mean[, j], lty = 2)
  points(parent_cutoffs, observed_prob[, j],
    pch = 16, cex = 0.6, col = "#B5261E"
  )
  if (j == 3) {
    legend(
      "topright",
      c(
        "Observed", "80% null interval",
        "95% null interval"
      ),
      col = c("#B5261E", "black", NA),
      fill = c(NA, "#777777", "#D0D0D0"), border = NA,
      lty = c(1, NA, NA), pch = c(16, NA, NA), bty = "n", cex = 0.9,
      seg.len = 1, x.intersp = 0.6, pt.cex = 0.7
    )
  }
  mtext(letters[j + 3], side = 3, adj = 0, line = 0.1, font = 2, cex = 1.25)
}
