## Recompute the course report's comparison table on exactly the same simulated data
## (same seed and simulation order as the report), with the original and the corrected estimators.
## Output: results/single_dataset.csv
source("R/load.R"); source("analysis/legacy_report_methods.R")
set.seed(12345)
sim <- lapply(REGIMES, function(th) do.call(simulate_jump_diffusion, c(list(n_obs = 2000), as.list(th))))

rows <- list()
for (reg in names(REGIMES)) {
  R <- sim[[reg]]$R
  set.seed(789); g_old <- gibbs_report(R)
  set.seed(2025); g_new <- fit_jd_gibbs(R)$summary
  st <- fit_jd_stan(R, chains = 4, iter = 2000, warmup = 1000, seed = 2025)
  est <- list(
    report_MLE   = fit_jd_mle_naive(R)$theta,
    report_EM    = fit_em_report(R),
    report_Gibbs = g_old,
    MLE          = fit_jd_mle(R)$theta,
    EM           = fit_jd_em(R)$theta,
    Gibbs_median = g_new[, "median"],
    Stan_median  = st$summary[, "median"])
  for (m in names(est)) rows[[length(rows) + 1]] <- data.frame(regime = reg, method = m, t(est[[m]][PARAM_NAMES]),
                                                              jumps_in_sample = sum(sim[[reg]]$Y))
  cat(sprintf("%s: Stan divergences %d, max R-hat %.3f; Gibbs max R-hat %.3f\n", reg, st$divergences, st$max_rhat,
              max(g_new[, "rhat"])))
}
out <- do.call(rbind, rows)
truth <- do.call(rbind, lapply(names(REGIMES), function(r) data.frame(regime = r, method = "TRUE", t(REGIMES[[r]]),
                                                                       jumps_in_sample = sum(sim[[r]]$Y))))
out <- rbind(truth, out); out <- out[order(match(out$regime, names(REGIMES))), ]
write.csv(out, "results/single_dataset.csv", row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)
