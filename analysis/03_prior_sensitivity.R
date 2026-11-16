## Prior sensitivity in the low regime (rare jumps): the Gibbs sampler on the report's low-regime data
## under three priors for the jump probability p, all other priors unchanged.
##   Beta(1, 9)  - default, prior mean 0.10
##   Beta(1, 49) - prior mean 0.02, close to the true value 0.01
##   Beta(1, 1)  - uniform
## Under the uniform prior the sampler occasionally wanders to large p (the 'fake jump' region), so every
## prior gets 4 long chains; short runs give unstable answers (check the rhat column).
## Usage: Rscript analysis/03_prior_sensitivity.R. Output: results/prior_sensitivity.csv
source("R/load.R")
set.seed(12345)
sim <- lapply(REGIMES, function(th) do.call(simulate_jump_diffusion, c(list(n_obs = 2000), as.list(th))))
R <- sim$low$R

priors <- list(`Beta(1, 9)` = c(1, 9), `Beta(1, 49)` = c(1, 49), `Beta(1, 1)` = c(1, 1))
rows <- list()
for (nm in names(priors)) {
  prior <- modifyList(GIBBS_PRIOR, list(a_p = priors[[nm]][1], b_p = priors[[nm]][2]))
  set.seed(2024)
  g <- fit_jd_gibbs(R, n_iter = 40000, burnin = 5000, thin = 10, chains = 4, prior = prior)$summary
  for (j in c("p", "mu_J", "sigma_J"))
    rows[[length(rows) + 1]] <- data.frame(prior_for_p = nm, parameter = j, truth = REGIMES$low[[j]],
      median = g[j, "median"], mean = g[j, "mean"], lo_2.5 = g[j, "lo"], hi_97.5 = g[j, "hi"], rhat = g[j, "rhat"])
}
out <- do.call(rbind, rows)
out$jumps_in_sample <- sum(sim$low$Y)
write.csv(out, "results/prior_sensitivity.csv", row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)
