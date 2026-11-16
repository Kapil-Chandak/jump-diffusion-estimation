## Bayesian estimation with Stan (Hamiltonian Monte Carlo / NUTS) on the marginalised likelihood.
## The model file is stan/jump_diffusion.stan; priors match the Gibbs sampler.

.stan_model_cache <- new.env()

jd_stan_model <- function(stan_file = file.path(.repo_root(), "stan", "jump_diffusion.stan")) {
  if (is.null(.stan_model_cache$m)) {
    rstan::rstan_options(auto_write = TRUE)
    .stan_model_cache$m <- rstan::stan_model(stan_file)
  }
  .stan_model_cache$m
}

#' Fit with Stan; returns posterior summaries and the main convergence diagnostics.
fit_jd_stan <- function(R, chains = 4, iter = 2000, warmup = 1000, seed = 2025, cores = 1) {
  fit <- rstan::sampling(jd_stan_model(), data = list(T = length(R), R = as.vector(R)), chains = chains,
                         iter = iter, warmup = warmup, seed = seed, cores = cores, refresh = 0)
  draws <- as.matrix(fit, pars = PARAM_NAMES)
  s <- rstan::summary(fit, pars = PARAM_NAMES)$summary
  summ <- cbind(mean = colMeans(draws), median = apply(draws, 2, median), sd = apply(draws, 2, sd),
                lo = apply(draws, 2, quantile, .025), hi = apply(draws, 2, quantile, .975), rhat = s[PARAM_NAMES, "Rhat"])
  list(summary = summ, draws = draws,
       divergences = sum(rstan::get_divergent_iterations(fit)), max_rhat = max(s[, "Rhat"], na.rm = TRUE))
}
