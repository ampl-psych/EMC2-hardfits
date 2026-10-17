## One page of small multiples: per case, the group mean whose chain means are furthest apart (in within-chain sds).
##   Rscript grid.R <dir with extracts> <out.png> <rows> <cols>
a <- commandArgs(TRUE); fs <- list.files(a[1], "[.]rds$", full.names = TRUE); nr <- as.integer(a[3]); nc <- as.integer(a[4])
png(a[2], 400 * nc, 340 * nr); par(mfrow = c(nr, nc), mar = c(2, 3, 3, 1))
for (f in fs) { x <- readRDS(f); ch <- x$chains; mu <- lapply(ch, `[[`, "theta_mu")
  sp <- apply(sapply(mu, rowMeans), 1, function(v) diff(range(v))) / apply(sapply(mu, function(m) apply(m, 1, sd)), 1, mean); j <- which.max(sp)
  sm <- sapply(ch, function(c) rowMeans(c$subj_ll) / attr(c$subj_ll, "scale")); ssp <- apply(sm, 1, function(v) diff(range(v)))
  matplot(ch[[1]]$idx, sapply(mu, function(m) m[j, ]), type = "l", lty = 1, col = seq_along(ch), xlab = "", ylab = "",
    main = sprintf("%s\n%s (spread %.1f sd) | subj spread max %.1f (%s)", x$meta$label, rownames(mu[[1]])[j], sp[j], max(ssp), names(ssp)[which.max(ssp)]), cex.main = 1.1) }
invisible(dev.off())
