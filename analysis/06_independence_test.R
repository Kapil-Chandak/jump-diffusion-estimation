## Does the independence test work when the model is true? Simulate many datasets per regime
## (jumps independent by construction), fit EM, flag days with P(jump | return) > 0.5 and test whether
## the flagged days are spread evenly over 8 blocks of 250 days. A valid test rejects in about 5% of
## datasets at the 5% level. The same test is applied to real data in analysis/07_real_data.R.
## Usage: Rscript analysis/06_independence_test.R [reps=200] [cores=2]
## Output: results/independence_test_simulated.csv
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) > 0) as.integer(args[1]) else 200
CORES <- if (length(args) > 1) as.integer(args[2]) else 2
source("R/load.R")
block <- rep(1:8, each = 250)

one <- function(rep) {
  set.seed(60600 + rep)                        # one seed per replicate: same results on any number of cores
  do.call(rbind, lapply(names(REGIMES), function(reg) {
    R <- do.call(simulate_jump_diffusion, c(list(n_obs = 2000), as.list(REGIMES[[reg]])))$R
    res <- jd_independence_test(fit_jd_em(R)$theta, R, block, B = 1999)
    data.frame(rep = rep, regime = reg, flagged = res[["flagged"]], p_value = res[["p_value"]])
  }))
}
raw <- do.call(rbind, par_map(seq_len(REPS), one, CORES))
out <- do.call(rbind, lapply(names(REGIMES), function(reg) {
  d <- raw[raw$regime == reg, ]
  data.frame(regime = reg, datasets = nrow(d), testable = sum(!is.na(d$p_value)),
             median_flagged_days = median(d$flagged), rejection_rate_5pct = mean(d$p_value < 0.05, na.rm = TRUE))
}))
write.csv(out, "results/independence_test_simulated.csv", row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
