# EMC2 hard fits: chain extracts

Compact extracts of hierarchical EMC2 fits that did and did not converge, for calibrating chain diagnostics
(chain culling, flat-region tests) offline. **No data are included**: each file holds sampler output only (draws of
group-level parameters, per-subject parameters and per-subject log-likelihoods), never the response times or choices.

`index.csv` lists every file: dataset, model, size, whether the fit met its stop rule, Rhat, the worst group-level
parameter, minutes, the sampler build, and a verdict read from the chain plots.

## Folders

| folder | what | build |
|---|---|---|
| `conflicteams/` | the 22 ConflictEAMs cases (ANT, Simon and Stroop conflict tasks) that never met the stop rule under the R sampler, fitted fresh on the C++ sampler, plus 3 easy models that converged (`Simon_LNR_MiSi`, `ANT_LNR_MiSi`, `Stroop_RDM_ViSi`) as controls. 6 of the 22 met the stop rule on this build and are controls too. | C++ sampler, EMC2 `d0b6cf0b` (`joint_move`) |
| `jm1/` | the same 22 ConflictEAMs cases on the next build, which adds a joint location move; each case was run with three seeds (`_c1`-`_c3`). 31 of the 60 runs are here: all three seeds of `ANT_LBA_BciVci`, `ANT_RDM_BcSi`, `ANT_RDM_Sci`, `ANT_RDM_VcSc`, `Simon_LBA_BciSci` and `Simon_WDM_BiVcSc`, the three DDMs (one seed each), and one seed (the one with the median group-mean Rhat, except `ANT_RDM_BiVi_c1`, the seed in which subject s543 sits in two modes) of `BcVc`, `BcVi`, `Bci`, `BiSi`, `BiVc`, `BiVi`, `VcSi`, `Vci`, `ViSi` and `Simon_LNR_MciSci`; `BcSc`, `BiSc` and `ViSc` are left out (they behave the same on both builds). 11 of the 31 met the stop rule (`jm1/stop_rule.csv`: outcome, checks, minutes and the final Rhats of each run). Three seeds show what one run cannot: a few ANT subjects (s543, s243, s304 and others) have two posterior modes 20-70 log-likelihood units apart, and which mode a run settles in is chance, under both builds. | C++ sampler, EMC2 `3c1bcaf1` (`joint_move_pr1`) |
| `eisenberg/` | the Eisenberg Simon data (518 subjects, RDM): five runs of the rewritten R sampler (`h6d_r1`-`r5`) and three of the legacy R sampler (`legacy_r1`-`r3`). The legacy runs look converged (Rhat at most 1.04, chains agree) but are wrong: a Stan reference puts subject parameters up to 46 posterior sd away. A diagnostic that only compares chains cannot see this. | R sampler (dev-sampler stage H6d), legacy R sampler |
| `plots/` | per file a four-panel chain plot (`plots/<folder>/<label>.png`) and one page per folder (`grid_*.png`) | |
| `tools/` | `extract.R` (writes the files), `diag_plots.R`, `grid.R`, `make_index.R` (builds `index.csv`) | |

The stop rule is `fit()`'s: in the sample stage, the largest group-mean Rhat and the .99 quantile of the subject
Rhats both below 1.1 (ConflictEAMs: up to 40 checks, 20 for DDMs). For Eisenberg, `converged` means the final
group-mean and subject Rhats are both below 1.1.

## File format

`readRDS("<folder>/<label>.rds")` gives `list(meta, chains)`.

`meta`: `label`, `dataset`, `build`, `model`, `type` (group-level structure, e.g. `"standard"`), `n_subjects`,
`n_pars`, `par_names`, `chain_n` (iterations per stage per chain in the saved fit), `rhat` (alpha Rhat per parameter,
maximum over subjects), `rhat_mu`, `rhat_Sigma` (as `gd_summary(..., stat_only = TRUE)` gave them).

`chains`: one element per chain, with

| element | content |
|---|---|
| `stage`, `idx` | stage and iteration number of each kept group-level draw (up to 1500, evenly spaced) |
| `theta_mu` | group means, parameters x draws |
| `theta_var` | group covariance, parameters x parameters x draws |
| `subj_ll_idx`, `subj_ll` | iteration numbers and per-subject log-likelihoods, subjects x draws, **stored as integers**: divide by `attr(subj_ll, "scale")` (1000) |
| `alpha_idx`, `alpha` | iteration numbers and subject parameters, parameters x subjects x draws, **stored as integers**: divide by `attr(alpha, "scale")` (10000) |
| `stage_all`, `idx_all` | stage and iteration number of every burn/adapt/sample iteration in the saved fit |
| `scale_move` | the sampler's joint scale/location move settings at the end of the fit (C++ build only) |

`subj_ll` and `alpha` are thinned more than the group level to keep each file near 15 MB (fits with 500+ subjects keep
about 500 log-likelihood draws and 90-140 alpha draws; the smaller Simon fits keep more). All three index vectors
refer to the same iteration numbering, so the components can be aligned. Parameters are on the sampled
(transformed) scale.

**The C++-build fits keep only their sample stage**: the saved objects hold no burn or adapt iterations, so
`stage` is `"sample"` throughout in `conflicteams/` and `jm1/`. The Eisenberg fits keep burn, adapt and sample.

## Verdicts

Read from the plots in `plots/` (chains in colour, chain 1 black, 2 red, 3 green, 4 blue, 5 cyan, 6 magenta):
"fine"; "slow drift" (all chains wander together over a group mean, no chain separate); "stuck chain" (one chain
held apart from the others for most of the run); "stuck subject(s)" (named subjects whose log-likelihood differs
between chains, i.e. separate modes).
