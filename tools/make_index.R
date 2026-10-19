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
## jm1/: three-seed runs on the next build; converged = fit()'s stop rule as the run ended (jm1/stop_rule.csv).
jm1_verdict <- c(
  ANT_LBA_BciVci_jm1_c1 = "slow drift (v_iMd; chains up to 2.9 within-chain sd apart)",
  ANT_LBA_BciVci_jm1_c2 = "slow drift (v_iMd)",
  ANT_LBA_BciVci_jm1_c3 = "stuck chains: chains 2 and 6 in a degenerate region from the burn stage on (sv_lMd -26 and -19, 2,700-4,600 log-lik units below the rest); the other four drift on v_iMd",
  ANT_RDM_BcSi_jm1_c1 = "stuck subject(s): s543 in the better mode in chains 2 and 4 only (66 log-lik units); group level fine",
  ANT_RDM_BcSi_jm1_c2 = "fine by the chains; every chain holds s543, s243 and s304 in the worse of their two modes (see ANT_RDM_BcSi_jm)",
  ANT_RDM_BcSi_jm1_c3 = "stuck subject(s): s543 and s243 in the better mode in chain 6 only (70 log-lik units); group level fine",
  ANT_RDM_BcVc_jm1_c2 = "slow drift (v; chains 4-5 high around iterations 1000-1500)",
  ANT_RDM_BcVi_jm1_c1 = "fine",
  ANT_RDM_Bci_jm1_c3 = "slow drift (v, v_iMd), mild",
  ANT_RDM_BiSi_jm1_c1 = "slow drift (v), mild",
  ANT_RDM_BiVc_jm1_c2 = "slow drift (v)",
  ANT_RDM_BiVi_jm1_c1 = "stuck chain: chain 2 offset on B_iMd for the whole run and in the worse mode for s543, s243, s454, s436 (20 log-lik units)",
  ANT_RDM_Sci_jm1_c1 = "stuck subject(s): s543, s243, s454, s436, s021 in two modes across chains (up to 61 log-lik units); plus slow drift (v; chains 3 within-chain sd apart)",
  ANT_RDM_Sci_jm1_c2 = "fine",
  ANT_RDM_Sci_jm1_c3 = "fine",
  ANT_RDM_VcSc_jm1_c1 = "stuck subject(s): s543, s021 and s254 in the better mode in chain 4 only (up to 73 log-lik units); plus slow drift (v)",
  ANT_RDM_VcSc_jm1_c2 = "fine",
  ANT_RDM_VcSc_jm1_c3 = "slow drift (v)",
  ANT_RDM_VcSi_jm1_c2 = "slow drift (v)",
  ANT_RDM_Vci_jm1_c3 = "slow drift (v_CIcong; chain 6 low over the last 1500 iterations)",
  ANT_RDM_ViSi_jm1_c3 = "slow drift (v)",
  Simon_DDM_Bi_jm1_c1 = "fine",
  Simon_DDM_BiVcSVc_jm1_c1 = "fine",
  Simon_LBA_BciSci_jm1_c1 = "slow drift along the B_CIRincong/B_iMd ridge at equal likelihood; chain 1 high on t0 over the last 1000 iterations",
  Simon_LBA_BciSci_jm1_c2 = "stuck chain: chain 3 offset on B_iMd/B_CIRincong at equal likelihood (group-level B/sv ridge)",
  Simon_LBA_BciSci_jm1_c3 = "slow drift along the B_CIRincong/B_iMd ridge at equal likelihood (chains 2.4 within-chain sd apart)",
  Simon_LNR_MciSci_jm1_c3 = "fine",
  Simon_WDM_BiVcSc_jm1_c1 = "fine at group level; subject Rhats on Z_lbCI up to 1.24 with no subject in a separate mode",
  Simon_WDM_BiVcSc_jm1_c2 = "fine at group level; subject Rhats on Z_lbCI up to 1.24 with no subject in a separate mode",
  Simon_WDM_BiVcSc_jm1_c3 = "fine",
  Stroop_DDM_BcSc_jm1_c1 = "slow drift (sv)")
if (file.exists("jm1/extract_log.csv")) {
  j <- read.csv("jm1/extract_log.csv", stringsAsFactors = FALSE)
  sr <- read.csv("jm1/stop_rule.csv", stringsAsFactors = FALSE)
  idx <- rbind(idx, data.frame(label = j$label, folder = "jm1", dataset = j$dataset, model = j$model, type = j$type,
    n_subjects = j$n_subjects, n_pars = j$n_pars, chains = j$chains, sample_iters = j$sample_iters,
    converged = ifelse(sr$met_stop_rule[match(j$label, sr$label)], "yes", "no"), max_rhat_alpha = j$max_rhat_alpha,
    max_rhat_mu = j$max_rhat_mu, max_rhat_sigma2 = j$max_rhat_Sigma, worst_parameter = j$worst_parameter,
    minutes = j$minutes, build = j$build, verdict = jm1_verdict[j$label], file_MB = j$file_MB, stringsAsFactors = FALSE))
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
