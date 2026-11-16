## Gibbs sampler for the jump-diffusion model, augmenting both the jump indicators Y_t and the
## jump sizes Z_t. With both augmented, every full conditional is a standard distribution:
##   (Y_t, Z_t) | rest : P(Y_t = 1) from the mixture; Z_t | Y_t = 1 normal (as in the EM E-step)
##   mu      | rest    : normal                (prior N(0, 0.05^2))
##   sigma^2 | rest    : inverse-gamma         (prior IG(2, 2e-4))
##   p       | rest    : Beta(a + n1, b + T - n1)   (prior Beta(1, 9))
##   mu_J    | rest    : normal                (prior N(0, 0.05^2)), using the n1 sampled jump sizes
##   sigma_J^2 | rest  : inverse-gamma         (prior IG(2, 2e-4))
## (The original report sampled an unconstrained two-component mixture and mapped it to these parameters
##  afterwards, so it did not sample from this model's posterior.)

GIBBS_PRIOR <- list(mu_sd = 0.05, s_shape = 2, s_rate = 2e-4, a_p = 1, b_p = 9, muJ_sd = 0.05, sJ_shape = 2, sJ_rate = 2e-4)

gibbs_jd_chain <- function(R, n_iter, burnin, thin, prior, init) {
  Tn <- length(R); mu <- init[["mu"]]; s2 <- init[["sigma"]]^2; p <- init[["p"]]
  mJ <- init[["mu_J"]]; sJ2 <- init[["sigma_J"]]^2
  keep <- seq(burnin + 1, n_iter, by = thin)
  out <- matrix(NA_real_, length(keep), 5, dimnames = list(NULL, PARAM_NAMES)); k <- 0
  for (it in seq_len(n_iter)) {
    tau <- jd_jump_prob(c(mu = mu, sigma = sqrt(s2), p = p, mu_J = mJ, sigma_J = sqrt(sJ2)), R)
    Y <- runif(Tn) < tau; n1 <- sum(Y); Z <- numeric(Tn)
    if (n1 > 0) { w <- sJ2 / (s2 + sJ2); Z[Y] <- rnorm(n1, mJ + w * (R[Y] - mu - mJ), sqrt(s2 * w)) }
    prec <- 1 / prior$mu_sd^2 + Tn / s2
    mu <- rnorm(1, (sum(R - Z) / s2) / prec, sqrt(1 / prec))
    s2 <- 1 / rgamma(1, prior$s_shape + Tn / 2, rate = prior$s_rate + 0.5 * sum((R - mu - Z)^2))
    p <- rbeta(1, prior$a_p + n1, prior$b_p + Tn - n1)
    precJ <- 1 / prior$muJ_sd^2 + n1 / sJ2
    mJ <- rnorm(1, (if (n1) sum(Z[Y]) else 0) / sJ2 / precJ, sqrt(1 / precJ))
    sJ2 <- 1 / rgamma(1, prior$sJ_shape + n1 / 2, rate = prior$sJ_rate + 0.5 * (if (n1) sum((Z[Y] - mJ)^2) else 0))
    if (it %in% keep) { k <- k + 1; out[k, ] <- c(mu, sqrt(s2), p, mJ, sqrt(sJ2)) }
  }
  out
}

#' Split-chain R-hat (Gelman-Rubin) for one parameter, from a list of chains.
rhat <- function(chains) {
  halves <- unlist(lapply(chains, function(x) { n <- floor(length(x) / 2); list(x[1:n], x[(n + 1):(2 * n)]) }), recursive = FALSE)
  n <- length(halves[[1]]); m <- length(halves)
  B <- n * var(sapply(halves, mean)); W <- mean(sapply(halves, var))
  sqrt(((n - 1) / n * W + B / n) / W)
}

#' Run several chains from dispersed starting values; return draws and posterior summaries.
fit_jd_gibbs <- function(R, n_iter = 6000, burnin = 2000, thin = 4, chains = 2, prior = GIBBS_PRIOR) {
  sr <- sd(R)
  inits <- list(c(mu = mean(R), sigma = sr / 1.5, p = 0.05, mu_J = -sr, sigma_J = 2 * sr),
                c(mu = mean(R), sigma = sr / 1.1, p = 0.01, mu_J = -3 * sr, sigma_J = sr))
  draws <- lapply(seq_len(chains), function(i) gibbs_jd_chain(R, n_iter, burnin, thin, prior, inits[[(i - 1) %% 2 + 1]]))
  all <- do.call(rbind, draws)
  summ <- t(apply(all, 2, function(v) c(mean = mean(v), median = median(v), sd = sd(v),
                                         lo = quantile(v, .025, names = FALSE), hi = quantile(v, .975, names = FALSE))))
  summ <- cbind(summ, rhat = sapply(PARAM_NAMES, function(j) rhat(lapply(draws, function(d) d[, j]))))
  list(summary = summ, draws = all)
}
