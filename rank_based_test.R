library(igraph)



# Script parameters
q <- round(seq(0.5, 0.9, by = 0.05), 2) # Sensitivity grid for Fisher's tests
q_plot <- c(0.5, 0.7, 0.9) # Figure and descriptive tree table
plot_indices <- match(q_plot, q)
seed <- 123
nrandom <- 10000
min_generations <- 2 # Maximum root-to-tip distance; index case is generation 0
min_size <- 20



# Load data
load("data/trees_processed.Rdata")
keep <- metadata$generations >= min_generations & metadata$size >= min_size
tree_indices <- which(keep)
metadata <- metadata[keep, , drop = FALSE]
trees <- trees[keep]
ntree <- length(trees)




# Function for creating a rooted tree from a Prufer sequence
make_rooted_from_prufer <- function(p, root, shuffle = FALSE) {
  if (shuffle) p <- sample(p)
  g <- make_from_prufer(p)
  parent <- bfs(g, root = root, parent = TRUE)$parent
  child <- setdiff(seq_len(vcount(g)), root)
  graph_from_edgelist(cbind(parent[child], child), directed = TRUE)
}



# Function for computing the rank-based statistic
transmission_rank_ratio <- function(g, q) {
  z <- degree(g, mode = "out")
  e <- ends(g, E(g), names = FALSE)

  r_child <- rank(z[e[, 2]], ties.method = "average") / nrow(e)
  r_parent <- rank(z[e[, 1]], ties.method = "average") / nrow(e)

  # Reuse the ranks across thresholds for the same tree.
  vapply(q, function(threshold) {
    top <- r_child >= threshold
    mean(r_parent[top]) / mean(r_parent[!top])
  }, numeric(1))
}



# Run the analysis
set.seed(seed)
nq <- length(q)

# Dimensions: observed [tree, q]; null [tree, q, replicate]
rho_obs <- matrix(NA, nrow = ntree, ncol = nq)
rho_null <- array(NA, dim = c(ntree, nq, nrandom))

for (i in seq_along(trees)) {
  message("Tree ", i, "/", length(trees), ": ", metadata$name[i])

  # Process the observed tree
  g <- upgrade_graph(trees[[i]]$tree)
  root <- which(degree(g, mode = "in") == 0)

  rho_obs[i, ] <- transmission_rank_ratio(g, q)

  # Evaluate each random tree at every q
  p <- to_prufer(as_undirected(g))

  for (k in seq_len(nrandom)) {
    g0 <- make_rooted_from_prufer(p, root, shuffle = TRUE)

    rho_null[i, , k] <- transmission_rank_ratio(g0, q)
  }
}



# Calculate null envelopes for the figure
plot_data <- lapply(plot_indices, function(j) {
  obs <- sort(rho_obs[, j])
  null <- rho_null[, j, ]
  null_sorted <- apply(null, 2, sort)

  envelope <- t(apply(null_sorted, 1, quantile,
    probs = c(0.025, 0.1, 0.5, 0.9, 0.975)
  ))

  list(obs = obs, envelope = envelope)
})



# Combine tree-level p-values using Fisher's statistic
tree_p <- matrix(NA,
  nrow = ntree,
  ncol = nq,
  dimnames = list(metadata$name, paste0("q_", q))
)
fisher <- matrix(NA, nrow = nq, ncol = nrandom + 1)

for (j in seq_along(q)) {
  ratio_all <- cbind(rho_obs[, j], rho_null[, j, ])

  p_all <- t(apply(ratio_all, 1, function(x) {
    rank(-x, ties.method = "max") / length(x)
  }))

  tree_p[, j] <- p_all[, 1]
  fisher[j, ] <- -2 * colSums(log(p_all))
}



# Calculate combined p-values from the joint randomizations
fisher_p <- t(apply(fisher, 1, function(x) {
  rank(-x, ties.method = "max") / length(x)
}))

global_test_table <- data.frame(
  test = paste0("q_", q),
  statistic = "fisher",
  observed_statistic = fisher[, 1],
  p_value = fisher_p[, 1],
  trees = ntree,
  randomizations = nrandom
)



# Create the tree table
tree_table <- data.frame(
  tree_index = tree_indices,
  tree_name = metadata$name,
  cases = metadata$size,
  disease = metadata$disease,
  country = metadata$country,
  year = metadata$year,
  citation = metadata$source,
  generations = metadata$generations,
  mean_offspring = vapply(trees, function(x) mean(x$out_degree), numeric(1)),
  mean_offspring_transmitting = sapply(trees, function(x) {
    z <- x$out_degree
    mean(z[z > 0])
  })
)

for (j in plot_indices) {
  # Tree-specific effect size and upper-tail Monte Carlo p-value, including
  # the observed configuration and counting ties in the upper tail.
  tree_table[[paste0("rho_q_", q[j])]] <- rho_obs[, j]
  tree_table[[paste0("p_upper_q_", q[j])]] <- tree_p[, j]

  ord <- order(rho_obs[, j])
  obs <- rho_obs[ord, j]
  null_sorted <- apply(rho_null[, j, ], 2, sort)

  # Across-tree order-statistic percentiles (not tree-specific p-values).
  # Calculate percentiles at observed ranks, with half weight for null ties
  null_percentile <- numeric(ntree)
  null_percentile[ord] <- rowMeans(null_sorted < obs) +
    0.5 * rowMeans(null_sorted == obs)

  # Average percentiles across tied observations
  null_percentile <- ave(null_percentile, rho_obs[, j], FUN = mean)
  tree_table[[paste0("null_percentile_q_", q[j])]] <- null_percentile
}

# Sort by decreasing null percentile
sort_q <- if (0.5 %in% q_plot) 0.5 else q_plot[1]
tree_table <- tree_table[
  order(tree_table[[paste0("null_percentile_q_", sort_q)]], decreasing = TRUE),
]
rownames(tree_table) <- NULL



# Save tables
dir.create("tables", showWarnings = FALSE)
write.csv(tree_table, "tables/trees.csv", row.names = FALSE)
write.csv(
  global_test_table, "tables/global_tests.csv", row.names = FALSE
)

# Export paste-ready Table S2 rows; the table and caption stay in ms.tex.
table_s2 <- tree_table[
  order(-tree_table$null_percentile_q_0.5, tree_table$disease),
]
sources <- read.csv("tables/tree_citations.csv", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(sources$tree_name))
source_keys <- sources$bibtex_keys[match(table_s2$tree_name, sources$tree_name)]
if (anyNA(source_keys) || any(!nzchar(source_keys))) {
  stop("Missing Table S2 citation keys in tables/tree_citations.csv")
}

latex_text <- function(x) {
  for (char in c("&", "%", "_", "#", "$")) {
    x <- gsub(char, paste0("\\", char), x, fixed = TRUE)
  }
  x
}
table_s2_rows <- paste0(
  latex_text(table_s2$disease), " & ",
  latex_text(paste(table_s2$country, table_s2$year, sep = ", ")), " & ",
  table_s2$cases, " & ",
  sprintf("%.2f", 100 * table_s2$null_percentile_q_0.5), " & ",
  sprintf("%.2f", 100 * table_s2$null_percentile_q_0.7), " & ",
  sprintf("%.2f", 100 * table_s2$null_percentile_q_0.9), " & ",
  "\\cite{", source_keys, "} \\\\"
)
table_s2_rows <- gsub("2005-07", "2005--07", table_s2_rows, fixed = TRUE)
writeLines(table_s2_rows, "tables/table_s2_rows.tex")



# Plot ranked observations and null intervals
xlim <- range(unlist(lapply(plot_data, function(x) c(x$obs, x$envelope))))

par(
  mar = c(3.8, 4.6, 1, 0.6),
  cex = 1,
  cex.axis = 1.1,
  cex.lab = 1.2,
  mgp = c(2, 0.7, 0)
)

for (j in seq_along(q_plot)) {
  obs <- plot_data[[j]]$obs
  envelope <- plot_data[[j]]$envelope
  position <- seq_along(obs)

  plot(NA,
    xlim = xlim,
    ylim = c(0, ntree),
    xlab = expression("Infector rank ratio, " * rho),
    ylab = "",
    bty = "l",
    yaxt = "n"
  )

  axis(2, las = 1)
  title(ylab = "Trees", line = 2.5)
  title(
    main = bquote(q == .(q_plot[j])),
    line = 0.1, cex.main = 1.1, font.main = 1
  )

  polygon(c(envelope[, 1], rev(envelope[, 5])),
    c(position, rev(position)),
    col = "#D0D0D0",
    border = NA
  )

  polygon(c(envelope[, 2], rev(envelope[, 4])),
    c(position, rev(position)),
    col = "#777777",
    border = NA
  )

  #lines(envelope[, 3], position, lty = 2)
  points(obs, position, pch = 16, cex = 0.6, col = "#B5261E")
  mtext(letters[j], side = 3, adj = 0, line = 0.1, font = 2, cex = 1.25)
}
