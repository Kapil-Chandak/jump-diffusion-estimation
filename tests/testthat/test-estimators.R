## Correctness checks for the likelihood and the estimators. Run: Rscript tests/run_tests.R
## (set RUN_STAN=1 to include the slower Stan check)

test_that("the mixture density integrates to one", {
  th <- REGIMES$medium
  f <- function(x) sapply(x, function(r) exp(jd_loglik(th, r)))
  expect_equal(integrate(f, -1, 1, subdivisions = 2000)$value, 1, tolerance = 1e-6)
})

test_that("EM never decreases the log-likelihood", {
  set.seed(1); th <- REGIMES$medium
  R <- do.call(simulate_jump_diffusion, c(list(n_obs = 2000), as.list(th)))$R
  ll <- sapply(1:60, function(k) fit_jd_em_single(R, maxit = k, tol = 0)$loglik)
  expect_true(all(diff(ll) > -1e-8))
})

test_that("EM and multi-start MLE reach the same maximum of the likelihood", {
  for (seed in 1:5) for (reg in names(REGIMES)) {
    set.seed(seed); th <- REGIMES[[reg]]
    R <- do.call(simulate_jump_diffusion, c(list(n_obs = 2000), as.list(th)))$R
    expect_lt(abs(fit_jd_em(R)$loglik - fit_jd_mle(R)$loglik), 0.01)
  }
})

test_that("EM recovers the parameters in a large sample", {
  set.seed(7); th <- REGIMES$medium
  R <- do.call(simulate_jump_diffusion, c(list(n_obs = 50000), as.list(th)))$R
  est <- fit_jd_em(R)$theta
  expect_lt(abs(est[["sigma"]] / th[["sigma"]] - 1), 0.03)
  for (j in c("p", "mu_J", "sigma_J")) expect_lt(abs(est[[j]] / th[[j]] - 1), 0.15)
})

test_that("Gibbs sampler converges and agrees with the EM estimate", {
  set.seed(11); th <- REGIMES$high
  R <- do.call(simulate_jump_diffusion, c(list(n_obs = 3000), as.list(th)))$R
  g <- fit_jd_gibbs(R, n_iter = 4000, burnin = 1000, thin = 2)
  em <- fit_jd_em(R)$theta
  expect_true(all(g$summary[, "rhat"] < 1.05))
  for (j in c("sigma", "p", "mu_J")) expect_lt(abs(g$summary[j, "median"] / em[[j]] - 1), 0.15)
})

test_that("Stan agrees with the Gibbs sampler (same model, same priors)", {
  skip_if(Sys.getenv("RUN_STAN") != "1", "set RUN_STAN=1 to run the Stan check")
  set.seed(11); th <- REGIMES$high
  R <- do.call(simulate_jump_diffusion, c(list(n_obs = 3000), as.list(th)))$R
  g <- fit_jd_gibbs(R, n_iter = 4000, burnin = 1000, thin = 2)$summary
  s <- fit_jd_stan(R, chains = 2)
  expect_equal(s$divergences, 0)
  for (j in c("sigma", "p", "mu_J", "sigma_J")) expect_lt(abs(s$summary[j, "median"] / g[j, "median"] - 1), 0.1)
})

test_that("the independence test accepts independent flags and rejects clustered ones", {
  set.seed(3); block <- rep(1:8, each = 250)
  pv <- replicate(200, jd_flag_test(runif(2000) < 0.03, block, B = 999)[["p_value"]])
  expect_lt(abs(mean(pv < 0.05) - 0.05), 0.04)                      # about 5% false rejections
  clustered <- block %in% c(3, 7) & runif(2000) < 0.1               # flags only in two "turbulent" blocks
  expect_lt(jd_flag_test(clustered, block, B = 999)[["p_value"]], 0.01)
})
