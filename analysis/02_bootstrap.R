## Bootstrap standard errors of the MLE and EM estimates, on the same simulated data as the report
## (seed 12345). In the low regime the SE of p is many times the estimate itself: the data barely
## identify the jump frequency.
## Usage: Rscript analysis/02_bootstrap.R [B=300]. Output: results/bootstrap_se.csv
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) > 0) as.integer(args[1]) else 300
source("R/load.R")
set.seed(12345)
sim <- lapply(REGIMES, function(th) do.call(simulate_jump_diffusion, c(list(n_obs = 2000), as.list(th))))

rows <- list()
for (reg in names(REGIMES)) {
  R <- sim[[reg]]$R
  for (m in c("MLE", "EM")) {
    fit_fun <- if (m == "MLE") fit_jd_mle else fit_jd_em
    est <- fit_fun(R)$theta
    set.seed(456)
    bs <- jd_bootstrap(R, fit_fun, B = B)
    rows[[length(rows) + 1]] <- data.frame(regime = reg, method = m, parameter = PARAM_NAMES,
      truth = REGIMES[[reg]][PARAM_NAMES], estimate = est[PARAM_NAMES], bootstrap_se = bs$se[PARAM_NAMES],
      se_over_estimate = abs(bs$se[PARAM_NAMES] / est[PARAM_NAMES]), failed_resamples = bs$failures, row.names = NULL)
  }
}
out <- do.call(rbind, rows)
write.csv(out, "results/bootstrap_se.csv", row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)
