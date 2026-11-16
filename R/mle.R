## Maximum likelihood by direct numerical optimisation.
## Parameters are optimised on an unconstrained scale (mu, log sigma, logit p, mu_J, log sigma_J).
## Two practical safeguards matter a great deal for this likelihood (see README, "Corrections"):
##   1. the data are rescaled to unit standard deviation, so the optimiser works on numbers of order one;
##   2. the optimiser is started from several points and the best maximum is kept, because the mixture
##      likelihood has several local maxima (including one where the "jump" component describes ordinary days).

.theta_from_phi <- function(phi, scale = 1) {
  c(mu = phi[1] * scale, sigma = exp(phi[2]) * scale, p = plogis(phi[3]),
    mu_J = phi[4] * scale, sigma_J = exp(phi[5]) * scale)
}

#' Single-start MLE exactly as in the original course report (kept to document the failure mode).
fit_jd_mle_naive <- function(R) {
  s0 <- sd(R) / 2
  phi0 <- c(mean(R), log(s0), qlogis(0.05), -0.03, log(s0))
  nll <- function(phi) { v <- -jd_loglik(.theta_from_phi(phi), R); if (is.finite(v)) v else 1e12 }
  o <- optim(phi0, nll, method = "BFGS")
  list(theta = .theta_from_phi(o$par), loglik = -o$value, converged = o$convergence == 0)
}

#' Recommended MLE: rescaled data, multiple starting values, best maximum kept.
fit_jd_mle <- function(R, p_starts = c(0.003, 0.01, 0.03, 0.1), jump_starts = c(-1, -3)) {
  sR <- sd(R); x <- R / sR
  nll <- function(phi) { v <- -jd_loglik(.theta_from_phi(phi), x); if (is.finite(v)) v else 1e12 }
  best <- NULL
  for (p0 in p_starts) for (m0 in jump_starts) {
    phi0 <- c(mean(x), log(0.8), qlogis(p0), m0, log(2))
    o <- tryCatch(optim(phi0, nll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  theta <- .theta_from_phi(best$par, sR)
  list(theta = theta, loglik = jd_loglik(theta, R), converged = best$convergence == 0)
}
