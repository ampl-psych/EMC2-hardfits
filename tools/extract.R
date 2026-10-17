## Write compact diagnostic extracts of saved EMC2 fits (one .rds per fit) for calibrating chain diagnostics offline.
##   Rscript extract.R <list.csv> <out dir> [cap MB = 15]
## list.csv columns: label, dataset, path (.RData or .rds holding one emc object), build, txt (optional fit record
## with minutes/status lines). Run with the EMC2 library that made the fits (EMC2_LIB).
## Per fit: burn/adapt/sample iterations only. theta_mu / theta_var at up to 1500 draws; subj_ll and alpha thinned
## separately so the file stays under the cap, each with its own iteration index; alpha stored as integers
## round(alpha * 1e4), subj_ll as round(subj_ll * 1e3) (attribute "scale"). Also writes <label>_trace.png (summed
## subject log-likelihood per chain) and a row of diagnostics used for the verdicts.
a <- commandArgs(TRUE); lst <- read.csv(a[1], stringsAsFactors = FALSE); out <- a[2]
cap <- if (length(a) > 2) as.numeric(a[3]) * 1e6 else 15e6
LIB <- Sys.getenv("EMC2_LIB"); EXTRA <- Sys.getenv("EXTRA_LIB")
.libPaths(c(if (nzchar(LIB)) LIB, if (nzchar(EXTRA)) EXTRA, .libPaths()))
suppressMessages(library(EMC2)); library(parallel)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

load_emc <- function(path) {
  if (grepl("\\.rds$", path)) return(readRDS(path))
  e <- new.env(); nm <- load(path, envir = e)
  for (n in nm) if (inherits(e[[n]], "emc")) return(e[[n]])
  e[[nm[1]]]
}
mx <- function(v) if (all(is.na(v))) NA else round(max(v, na.rm = TRUE), 3)
thin_idx <- function(keep, n) if (length(keep) > n) keep[round(seq(1, length(keep), length.out = n))] else keep
as_scaled <- function(x, s) { y <- array(as.integer(round(x * s)), dim(x), dimnames(x)); attr(y, "scale") <- s; y }

extract_one <- function(emc, label, budget_ll, budget_alpha) {
  s1 <- emc[[1]]$samples; n <- s1$idx
  keep <- which(s1$stage[seq_len(n)] %in% c("burn", "adapt", "sample"))
  nch <- length(emc); nsub <- dim(s1$alpha)[2]; npar <- dim(s1$alpha)[1]
  i_g <- thin_idx(keep, 1500)
  i_ll <- thin_idx(keep, max(50, min(1500, floor(budget_ll / (nch * nsub * 4)))))
  i_a <- thin_idx(keep, max(30, min(1500, floor(budget_alpha / (nch * nsub * npar * 4)))))
  ch <- lapply(emc, function(x) {
    s <- x$samples
    list(stage = s$stage[i_g], idx = i_g,
         theta_mu = s$theta_mu[, i_g, drop = FALSE],
         theta_var = if (!is.null(s$theta_var)) s$theta_var[, , i_g, drop = FALSE],
         subj_ll_idx = i_ll, subj_ll = as_scaled(s$subj_ll[, i_ll, drop = FALSE], 1e3),
         alpha_idx = i_a, alpha = as_scaled(s$alpha[, , i_a, drop = FALSE], 1e4),
         stage_all = s$stage[keep], idx_all = keep,
         scale_move = attr(s, "scale_move"), pm_settings = NULL)
  })
  ch
}

diagnose <- function(emc) {      # where the chains disagree, from the summed and per-subject log-likelihoods
  s1 <- emc[[1]]$samples; n <- s1$idx
  smp <- which(s1$stage[seq_len(n)] == "sample"); if (length(smp) < 20) smp <- which(s1$stage[seq_len(n)] == "adapt")
  smp <- smp[smp > smp[1] + length(smp) / 3]                       # last two thirds of the sample stage
  ll <- sapply(emc, function(x) colSums(x$samples$subj_ll[, smp, drop = FALSE]))     # draws x chains
  cm <- colMeans(ll); w <- mean(apply(ll, 2, sd))
  z_chain <- (cm - median(cm)) / w
  sm <- sapply(emc, function(x) rowMeans(x$samples$subj_ll[, smp, drop = FALSE]))     # subjects x chains
  sw <- sapply(emc, function(x) apply(x$samples$subj_ll[, smp, drop = FALSE], 1, sd))
  spread <- apply(sm, 1, function(v) max(v) - min(v)); ratio <- spread / rowMeans(sw)
  stuck <- which(spread > 10 & ratio > 4)
  list(z_chain = round(z_chain, 1), worst_chain_z = round(max(abs(z_chain)), 1),
       stuck_subjects = rownames(emc[[1]]$samples$subj_ll)[stuck][order(-spread[stuck])],
       max_subject_spread = round(max(spread), 1), ll = ll, smp = smp)
}

run <- function(r) {
  t0 <- Sys.time()
  emc <- tryCatch(load_emc(r$path), error = function(e) e)
  if (inherits(emc, "error")) return(data.frame(label = r$label, error = conditionMessage(emc)))
  rh_a <- tryCatch(unlist(gd_summary(emc, selection = "alpha", stat_only = TRUE)), error = function(e) NA)
  rh_m <- tryCatch(unlist(gd_summary(emc, selection = "mu", stat_only = TRUE)), error = function(e) NA)
  rh_s <- tryCatch(unlist(gd_summary(emc, selection = "sigma2", stat_only = TRUE)), error = function(e) NA)
  mdl <- tryCatch(emc[[1]]$model()$type, error = function(e) tryCatch(emc[[1]]$model$type, error = function(e) NA))
  meta <- list(label = r$label, dataset = r$dataset, build = r$build, model = mdl,
               type = emc[[1]]$type, n_subjects = dim(emc[[1]]$samples$alpha)[2], n_pars = length(emc[[1]]$par_names),
               par_names = emc[[1]]$par_names, chain_n = chain_n(emc), rhat = rh_a, rhat_mu = rh_m, rhat_Sigma = rh_s)
  bl <- 16e6; ba <- 40e6
  for (k in 1:4) {
    f <- file.path(out, paste0(r$label, ".rds"))
    ch <- extract_one(emc, r$label, bl, ba)
    saveRDS(list(meta = meta, chains = ch), f, compress = "xz")
    sz <- file.size(f); if (sz <= cap) break
    ba <- ba * cap / sz * .85; bl <- bl * min(1, cap / sz * .95)
  }
  d <- diagnose(emc)
  png(file.path(out, paste0(r$label, "_trace.png")), 1000, 450)
  matplot(d$smp, d$ll, type = "l", lty = 1, xlab = "iteration", ylab = "summed subject log-likelihood",
          main = paste(r$label, "- sample stage, last two thirds, by chain")); dev.off()
  txt <- if (!is.na(r$txt) && nzchar(r$txt) && file.exists(as.character(r$txt))) readLines(r$txt) else character(0)
  g <- function(k) { v <- sub(paste0("^", k, ": "), "", grep(paste0("^", k, ": "), txt, value = TRUE)); if (length(v)) v[length(v)] else NA }
  worst <- tryCatch({ v <- c(rh_m, setNames(rh_s, paste0("sigma2:", names(rh_s)))); names(v)[which.max(v)] }, error = function(e) NA)
  data.frame(label = r$label, dataset = r$dataset, model = meta$model, type = meta$type,
             n_subjects = meta$n_subjects, n_pars = meta$n_pars, chains = length(emc),
             sample_iters = chain_n(emc)[1, "sample"],
             max_rhat_alpha = mx(rh_a), q99_rhat_alpha = round(unname(quantile(rh_a, .99, na.rm = TRUE)), 3),
             max_rhat_mu = mx(rh_m), max_rhat_Sigma = mx(rh_s), worst_parameter = worst,
             minutes = as.numeric(g("minutes")), status = g("status"), build = r$build,
             worst_chain_z = d$worst_chain_z, chain_z = paste(d$z_chain, collapse = " "),
             max_subject_spread = d$max_subject_spread, n_stuck_subjects = length(d$stuck_subjects),
             stuck_subjects = paste(head(d$stuck_subjects, 8), collapse = " "),
             file_MB = round(file.size(file.path(out, paste0(r$label, ".rds"))) / 1e6, 1),
             group_draws = length(ch[[1]]$idx), subj_ll_draws = length(ch[[1]]$subj_ll_idx),
             alpha_draws = length(ch[[1]]$alpha_idx), extract_secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
rows <- mclapply(split(lst, seq_len(nrow(lst))), function(r) tryCatch(run(r), error = function(e) data.frame(label = r$label, error = conditionMessage(e))),
                 mc.cores = as.integer(Sys.getenv("MC_CORES", Sys.getenv("PBS_NCPUS", "2"))), mc.preschedule = FALSE)
cols <- unique(unlist(lapply(rows, names)))
res <- do.call(rbind, lapply(rows, function(x) { for (k in setdiff(cols, names(x))) x[[k]] <- NA; x[cols] }))
write.csv(res, file.path(out, "extract_log.csv"), row.names = FALSE)
print(res)
