## index.csv from the extract logs plus the stop-rule outcome and the plot verdicts (chains numbered 1..n as stored;
## in the plots chain 1 black, 2 red, 3 green, 4 blue, 5 cyan, 6 magenta).
conv <- c("ANT_RDM_BcSi_jm", "ANT_RDM_ViSc_jm", "ANT_RDM_BiSc_jm", "Simon_WDM_BiVcSc_jm", "Stroop_DDM_BcSc_jm",
          "Simon_DDM_Bi_jm", "Simon_LNR_MiSi_jm", "ANT_LNR_MiSi_jm", "Stroop_RDM_ViSi_jm")
verdict <- c(
  ANT_LBA_BciVci_jm = "slow drift (v_iMd)",
  ANT_LNR_MiSi_jm = "fine (control)",
  ANT_RDM_BcSc_jm = "slow drift (v; early transient in chains 4-5)",
  ANT_RDM_BcSi_jm = "fine",
  ANT_RDM_BcVc_jm = "slow drift (v)",
  ANT_RDM_BcVi_jm = "slow drift (v)",
  ANT_RDM_Bci_jm = "slow drift (v), mild",
  ANT_RDM_BiSc_jm = "fine",
  ANT_RDM_BiSi_jm = "slow drift (v; chains separated by 2.5 within-chain sd)",
  ANT_RDM_BiVc_jm = "slow drift (v_CIcong; chain 4 low for ~1500 iterations)",
  ANT_RDM_BiVi_jm = "slow drift (v), mild",
  ANT_RDM_Sci_jm = "stuck subject(s): s304, s254, s243, s436 hop between two modes (up to 14 log-lik units apart); plus slow drift (v)",
  ANT_RDM_VcSc_jm = "slow drift (v)",
  ANT_RDM_VcSi_jm = "stuck chain: chain 1 in the better mode for s254, s543, s454, s304 (27-50 log-lik units), chains 2-6 in the worse one; plus slow drift (v)",
  ANT_RDM_Vci_jm = "stuck chain: chain 1 offset on v_iMd for most of the run; slow drift",
  ANT_RDM_ViSc_jm = "fine",
  ANT_RDM_ViSi_jm = "slow drift (v)",
  Simon_DDM_Bi_jm = "fine",
  Simon_DDM_BiVcSVc_jm = "slow drift (sv_CIcong)",
  Simon_LBA_BciSci_jm = "stuck chain: chain 2 offset on B_CIRincong at equal likelihood (group-level B/sv ridge)",
  Simon_LNR_MciSci_jm = "stuck chain / slow drift: chain 2 high on m_iMd for most of the run",
  Simon_LNR_MiSi_jm = "fine (control)",
  Simon_WDM_BiVcSc_jm = "fine after a transient: the sample stage starts with Z_lbCI still falling (first ~400 iterations)",
  Stroop_DDM_BcSc_jm = "fine",
  Stroop_RDM_ViSi_jm = "fine (control)")
x <- read.csv("conflicteams/extract_log.csv", stringsAsFactors = FALSE)
idx <- data.frame(label = x$label, folder = "conflicteams", dataset = x$dataset, model = x$model, type = x$type,
  n_subjects = x$n_subjects, n_pars = x$n_pars, chains = x$chains, sample_iters = x$sample_iters,
  converged = ifelse(x$label %in% conv, "yes", "no"), max_rhat_alpha = x$max_rhat_alpha, max_rhat_mu = x$max_rhat_mu,
  max_rhat_sigma2 = x$max_rhat_Sigma, worst_parameter = x$worst_parameter, minutes = x$minutes, build = x$build,
  verdict = verdict[x$label], file_MB = x$file_MB, stringsAsFactors = FALSE)
if (file.exists("eisenberg/extract_log.csv")) {
  e <- read.csv("eisenberg/extract_log.csv", stringsAsFactors = FALSE)
  ev <- if (file.exists("tools/eisenberg_verdicts.csv")) read.csv("tools/eisenberg_verdicts.csv", stringsAsFactors = FALSE) else data.frame(label = e$label, converged = NA, verdict = NA)
  idx <- rbind(idx, data.frame(label = e$label, folder = "eisenberg", dataset = e$dataset, model = e$model, type = e$type,
    n_subjects = e$n_subjects, n_pars = e$n_pars, chains = e$chains, sample_iters = e$sample_iters,
    converged = ev$converged[match(e$label, ev$label)], max_rhat_alpha = e$max_rhat_alpha, max_rhat_mu = e$max_rhat_mu,
    max_rhat_sigma2 = e$max_rhat_Sigma, worst_parameter = e$worst_parameter, minutes = e$minutes, build = e$build,
    verdict = ev$verdict[match(e$label, ev$label)], file_MB = e$file_MB, stringsAsFactors = FALSE))
}
## worst group-level parameter by split-Rhat on the extract's sample-stage draws (group means mu: and log group variances lsd:)
split_rhat <- function(m) {                    # m: draws x chains
  h <- floor(nrow(m) / 2); m <- cbind(m[1:h, , drop = FALSE], m[(h + 1):(2 * h), , drop = FALSE])
  W <- mean(apply(m, 2, var)); B <- h * var(colMeans(m)); sqrt(((h - 1) / h * W + B / h) / W)
}
worst <- function(f) {
  ch <- readRDS(f)$chains; k <- which(ch[[1]]$stage == "sample")
  mu <- sapply(rownames(ch[[1]]$theta_mu), function(p) split_rhat(sapply(ch, function(c) c$theta_mu[p, k])))
  lv <- if (!is.null(ch[[1]]$theta_var)) sapply(rownames(ch[[1]]$theta_mu), function(p)
    split_rhat(sapply(ch, function(c) log(c$theta_var[p, p, k])))) else NULL
  v <- c(setNames(mu, paste0("mu:", names(mu))), if (!is.null(lv)) setNames(lv, paste0("lsd:", names(lv))))
  sprintf("%s (%.2f)", names(v)[which.max(v)], max(v))
}
idx$worst_parameter <- vapply(file.path(idx$folder, paste0(idx$label, ".rds")), worst, "")
write.csv(idx, "index.csv", row.names = FALSE)
