dir.create("figures", showWarnings = FALSE)
pdf("figures/empirical_combined.pdf", width = 8, height = 6, useDingbats = FALSE)
tryCatch({
  par(mfrow = c(2, 3))
  source("rank_based_test_v2.R", local = new.env())
  source("dyad_analysis.R", local = new.env())
}, finally = dev.off())
