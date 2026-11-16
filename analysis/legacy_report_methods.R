## The original course-report estimators, kept unchanged (apart from formatting) so that the
## corrections can be demonstrated on the same data. Do not use these for estimation.
##  - fit_em_report: EM for an UNCONSTRAINED two-component normal mixture, then
##    sigma_J^2 = max(v1 - v0, 0). Not the maximum likelihood estimate of the jump-diffusion model.
##  - gibbs_report: Gibbs sampler for the same unconstrained mixture, mapped afterwards.

fit_em_report <- function(R, max_iter = 500, tol = 1e-6) {
  Tlen <- length(R); R_mean <- mean(R); R_sd <- sd(R)
  pi <- 0.05; m0 <- R_mean; m1 <- R_mean + sign(-R_mean) * R_sd; if (!is.finite(m1)) m1 <- R_mean - R_sd
  v0 <- max(R_sd^2 / 2, 1e-6); v1 <- max(R_sd^2, 1e-6); ll_old <- -Inf
  for (iter in seq_len(max_iter)) {
    f0 <- dnorm(R, m0, sqrt(v0)); f1 <- dnorm(R, m1, sqrt(v1))
    g <- pi * f1 / (pi * f1 + (1 - pi) * f0); g[!is.finite(g)] <- 0.5
    N1 <- sum(g); N0 <- Tlen - N1; pi_n <- N1 / Tlen
    m0_n <- if (N0 > 0) sum((1 - g) * R) / N0 else m0; m1_n <- if (N1 > 0) sum(g * R) / N1 else m1
    v0_n <- max(if (N0 > 0) sum((1 - g) * (R - m0_n)^2) / N0 else v0, 1e-8)
    v1_n <- max(if (N1 > 0) sum(g * (R - m1_n)^2) / N1 else v1, 1e-8)
    ll <- sum(log((1 - pi_n) * dnorm(R, m0_n, sqrt(v0_n)) + pi_n * dnorm(R, m1_n, sqrt(v1_n))))
    done <- abs(ll - ll_old) < tol * (1 + abs(ll_old))
    pi <- pi_n; m0 <- m0_n; m1 <- m1_n; v0 <- v0_n; v1 <- v1_n; ll_old <- ll
    if (done) break
  }
  c(mu = m0, sigma = sqrt(max(v0, 1e-8)), p = pi, mu_J = m1 - m0, sigma_J = sqrt(max(v1 - v0, 0)))
}

gibbs_report <- function(R, n_iter = 10000, burnin = 2000, thin = 10, a_p = 1, b_p = 9,
                         m0_pr = 0, kappa0 = 0.01, alpha0 = 2, beta0 = 2e-4) {
  Tlen <- length(R); keep <- seq(burnin + 1, n_iter, by = thin)
  draws <- matrix(NA_real_, length(keep), 5, dimnames = list(NULL, PARAM_NAMES)); k <- 0
  p <- 0.05; m0 <- mean(R); m1 <- mean(R) - sd(R); v0 <- max(sd(R)^2 / 2, 1e-6); v1 <- max(sd(R)^2, 1e-6)
  Z <- rbinom(Tlen, 1, p)
  upd <- function(x, n) {
    if (n == 0) { v <- 1 / rgamma(1, alpha0, rate = beta0); return(c(rnorm(1, m0_pr, sqrt(v / kappa0)), v)) }
    xb <- mean(x); kn <- kappa0 + n; mn <- (kappa0 * m0_pr + n * xb) / kn
    bn <- beta0 + 0.5 * sum((x - xb)^2) + kappa0 * n * (xb - m0_pr)^2 / (2 * kn)
    v <- max(1 / rgamma(1, alpha0 + n / 2, rate = bn), 1e-8); c(rnorm(1, mn, sqrt(v / kn)), v)
  }
  for (it in seq_len(n_iter)) {
    f0 <- dnorm(R, m0, sqrt(v0)); f1 <- dnorm(R, m1, sqrt(v1))
    pr <- p * f1 / (p * f1 + (1 - p) * f0); pr[!is.finite(pr)] <- 0.5
    Z <- rbinom(Tlen, 1, pr); n1 <- sum(Z); n0 <- Tlen - n1
    p <- rbeta(1, a_p + n1, b_p + n0)
    a <- upd(R[Z == 0], n0); m0 <- a[1]; v0 <- a[2]
    b <- upd(R[Z == 1], n1); m1 <- b[1]; v1 <- b[2]
    if (it %in% keep) { k <- k + 1; draws[k, ] <- c(m0, sqrt(v0), p, m1 - m0, sqrt(max(v1 - v0, 0))) }
  }
  colMeans(draws)
}
