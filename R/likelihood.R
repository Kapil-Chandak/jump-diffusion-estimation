## Log-likelihood of the jump-diffusion model.
## Marginally, R_t is a two-component normal mixture:
##   (1 - p) N(mu, sigma^2) + p N(mu + mu_J, sigma^2 + sigma_J^2).
## Note the constraint: the jump component's variance is the diffusion variance PLUS sigma_J^2.

#' Per-observation log-densities of the two components (with mixing weights included).
jd_log_components <- function(theta, R) {
  s2 <- theta[["sigma"]]^2
  cbind(no_jump = log1p(-theta[["p"]]) + dnorm(R, theta[["mu"]], sqrt(s2), log = TRUE),
        jump    = log(theta[["p"]]) + dnorm(R, theta[["mu"]] + theta[["mu_J"]], sqrt(s2 + theta[["sigma_J"]]^2), log = TRUE))
}

#' Total log-likelihood, evaluated with log-sum-exp for numerical stability.
jd_loglik <- function(theta, R) {
  lc <- jd_log_components(theta, R)
  m <- pmax(lc[, 1], lc[, 2])
  sum(m + log(exp(lc[, 1] - m) + exp(lc[, 2] - m)))
}

#' Posterior probability that each day is a jump day.
jd_jump_prob <- function(theta, R) {
  lc <- jd_log_components(theta, R)
  1 / (1 + exp(lc[, 1] - lc[, 2]))
}
