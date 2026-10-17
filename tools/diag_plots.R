## Four-panel chain plots from the extracts (no EMC2 needed):   Rscript diag_plots.R <dir with extracts> <png dir>
## 1 summed subject log-likelihood by chain; 2 the group mean with the largest between-chain spread by chain;
## 3 per subject, the spread of the chain means of the log-likelihood (sorted, top 30); 4 that worst subject's
## log-likelihood by chain. Draws are the extract's (thinned) iterations; chains in colour.
a <- commandArgs(TRUE); dir.create(a[2], showWarnings = FALSE)
for (f in list.files(a[1], "\\.rds$", full.names = TRUE)) {
  x <- readRDS(f); ch <- x$chains; lab <- x$meta$label; nc <- length(ch)
  ll <- sapply(ch, function(c) colSums(c$subj_ll) / attr(c$subj_ll, "scale"))
  it <- ch[[1]]$subj_ll_idx
  mu <- lapply(ch, function(c) c$theta_mu)
  spread_mu <- apply(sapply(mu, rowMeans), 1, function(v) diff(range(v))) / apply(sapply(mu, function(m) apply(m, 1, sd)), 1, mean)
  j <- which.max(spread_mu)
  sm <- sapply(ch, function(c) rowMeans(c$subj_ll) / attr(c$subj_ll, "scale"))
  sp <- apply(sm, 1, function(v) diff(range(v))); o <- order(-sp); w <- o[1]
  png(file.path(a[2], paste0(lab, ".png")), 1400, 900)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  matplot(it, ll, type = "l", lty = 1, col = seq_len(nc), xlab = "iteration", ylab = "summed subject log-lik",
          main = paste(lab, "| Rhat alpha max", round(max(x$meta$rhat, na.rm = TRUE), 2), "mu max", round(max(x$meta$rhat_mu, na.rm = TRUE), 2)))
  matplot(ch[[1]]$idx, sapply(mu, function(m) m[j, ]), type = "l", lty = 1, col = seq_len(nc), xlab = "iteration",
          ylab = rownames(mu[[1]])[j], main = paste("group mean with largest chain spread:", rownames(mu[[1]])[j]))
  barplot(sp[o[1:min(30, length(o))]], names.arg = names(sp)[o[1:min(30, length(o))]], las = 2, cex.names = .8,
          main = "between-chain spread of subject mean log-lik (top 30)", ylab = "log-lik units")
  matplot(it, sapply(ch, function(c) c$subj_ll[w, ] / attr(c$subj_ll, "scale")), type = "l", lty = 1, col = seq_len(nc),
          xlab = "iteration", ylab = "log-lik", main = paste("worst subject", names(sp)[w]))
  dev.off()
}
