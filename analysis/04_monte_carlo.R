## Monte Carlo study: the three regimes, T = 2000 days, many simulated datasets.
## Every replicate: original report MLE and EM (to document the corrections), corrected MLE and EM, Gibbs.
## Stan (slowest) on the first STAN_REPS replicates. Replicates are saved as they finish (resumable).
## Usage: Rscript analysis/04_monte_carlo.R [reps=200] [stan_reps=25] [cores=2]
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) > 0) as.integer(args[1]) else 200
STAN_REPS <- if (length(args) > 1) as.integer(args[2]) else 25
CORES <- if (length(args) > 2) as.integer(args[3]) else 2
source("R/load.R"); source("analysis/legacy_report_methods.R"); library(parallel)   # parallel: RNG streams
invisible(jd_stan_model())                            # compile once before forking
dir.create("results/mc_reps", showWarnings = FALSE, recursive = TRUE)
RNGkind("L'Ecuyer-CMRG"); set.seed(20260101)
seeds <- local({ s <- .Random.seed; lapply(seq_len(REPS), function(i) { s <<- nextRNGStream(s); s }) })

one_rep <- function(rep) {
  f <- sprintf("results/mc_reps/rep_%04d.rds", rep)
  if (file.exists(f)) return(readRDS(f))
  assign(".Random.seed", seeds[[rep]], envir = .GlobalEnv)
  rows <- list()
  add <- function(reg, method, th, lo = NULL, hi = NULL, njump, ll = NA)
    rows[[length(rows) + 1]] <<- data.frame(rep, regime = reg, method, parameter = PARAM_NAMES, estimate = unname(th[PARAM_NAMES]),
      truth = unname(REGIMES[[reg]][PARAM_NAMES]), lo = if (is.null(lo)) NA else unname(lo[PARAM_NAMES]),
      hi = if (is.null(hi)) NA else unname(hi[PARAM_NAMES]), jumps_in_sample = njump, loglik = ll)
  for (reg in names(REGIMES)) {
    sim <- do.call(simulate_jump_diffusion, c(list(n_obs = 2000), as.list(REGIMES[[reg]]))); R <- sim$R; nj <- sum(sim$Y)
    add(reg, "report_MLE", fit_jd_mle_naive(R)$theta, njump = nj)
    add(reg, "report_EM", fit_em_report(R), njump = nj)
    f_mle <- fit_jd_mle(R); f_em <- fit_jd_em(R)
    add(reg, "MLE", f_mle$theta, njump = nj, ll = f_mle$loglik)
    add(reg, "EM", f_em$theta, njump = nj, ll = f_em$loglik)
    g <- fit_jd_gibbs(R)$summary
    add(reg, "Gibbs", g[, "median"], g[, "lo"], g[, "hi"], nj)
    if (rep <= STAN_REPS) {
      s <- fit_jd_stan(R, chains = 2, iter = 2000, warmup = 1000, seed = rep)
      add(reg, "Stan", s$summary[, "median"], s$summary[, "lo"], s$summary[, "hi"], nj)
      rows[[length(rows)]]$divergences <- s$divergences
    }
  }
  out <- do.call(rbind, lapply(rows, function(r) { if (is.null(r$divergences)) r$divergences <- NA; r }))
  saveRDS(out, f); out
}
res <- do.call(rbind, par_map(seq_len(REPS), one_rep, CORES))
saveRDS(res, "results/monte_carlo.rds")
cat("done:", REPS, "replicates\n")
