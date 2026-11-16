## Exact EM algorithm for the jump-diffusion model.
## Complete data: returns R_t, jump indicators Y_t and jump sizes Z_t. Given Y_t = 1,
##   R_t - mu = Z_t + sigma * eps_t,  so  Z_t | R_t, Y_t = 1  ~  N(m_t, v),
##   m_t = mu_J + sigma_J^2 / (sigma^2 + sigma_J^2) * (R_t - mu - mu_J),   v = sigma^2 sigma_J^2 / (sigma^2 + sigma_J^2).
## Every M-step update is then in closed form, and the log-likelihood never decreases.
## (The original report instead fitted an unconstrained two-component mixture and derived sigma_J^2 as a
##  difference of variances; that is not the maximum likelihood estimate of this model.)

fit_jd_em_single <- function(R, p0 = 0.05, jump0 = -1, maxit = 5000, tol = 1e-10) {
  sr <- sd(R)
  mu <- mean(R); s <- sr / 1.2; p <- p0; mJ <- jump0 * sr; sJ <- 2 * sr
  ll_old <- -Inf; converged <- FALSE
  for (it in seq_len(maxit)) {
    th <- c(mu = mu, sigma = s, p = p, mu_J = mJ, sigma_J = sJ)
    tau <- jd_jump_prob(th, R)                                   # E-step: P(Y_t = 1 | R_t)
    w <- sJ^2 / (s^2 + sJ^2)
    m <- mJ + w * (R - mu - mJ); v <- s^2 * w                   # E-step: Z_t | R_t, Y_t = 1
    EYZ <- tau * m; EYZ2 <- tau * (m^2 + v)
    p  <- mean(tau)                                              # M-step
    mu <- mean(R - EYZ)
    s  <- sqrt(mean((R - mu)^2 - 2 * (R - mu) * EYZ + EYZ2))
    mJ <- sum(EYZ) / sum(tau)
    sJ <- sqrt(max(sum(EYZ2 - 2 * mJ * EYZ + tau * mJ^2) / sum(tau), 1e-12))
    if (p < 1e-5 || sJ < 1e-7) break                             # jump component vanished
    ll <- jd_loglik(c(mu = mu, sigma = s, p = p, mu_J = mJ, sigma_J = sJ), R)
    if (abs(ll - ll_old) < tol * (1 + abs(ll))) { converged <- TRUE; break }
    ll_old <- ll
  }
  theta <- c(mu = mu, sigma = s, p = p, mu_J = mJ, sigma_J = sJ)
  list(theta = theta, loglik = jd_loglik(theta, R), iterations = it, converged = converged)
}

#' EM from several starting values; the fit with the highest log-likelihood is kept.
fit_jd_em <- function(R, p_starts = c(0.003, 0.01, 0.03, 0.1), jump_starts = c(-1, -3), ...) {
  fits <- list()
  for (p0 in p_starts) for (j0 in jump_starts)
    fits[[length(fits) + 1]] <- tryCatch(fit_jd_em_single(R, p0, j0, ...), error = function(e) NULL)
  fits <- Filter(function(f) !is.null(f) && is.finite(f$loglik), fits)
  fits[[which.max(sapply(fits, `[[`, "loglik"))]]
}
