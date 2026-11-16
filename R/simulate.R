## Simulation of the discrete-time jump-diffusion model
##   R_t = mu + sigma * eps_t + Y_t * Z_t,   eps_t ~ N(0,1),  Y_t ~ Bernoulli(p),  Z_t ~ N(mu_J, sigma_J^2)

#' Simulate T daily log-returns from the jump-diffusion model.
#' Returns the returns together with the latent jump indicators and jump sizes, so that
#' estimates can be compared with the truth.
simulate_jump_diffusion <- function(n_obs, mu, sigma, p, mu_J, sigma_J) {
  Y   <- rbinom(n_obs, size = 1, prob = p)       # jump indicators
  Z   <- rnorm(n_obs, mean = mu_J, sd = sigma_J) # jump sizes
  eps <- rnorm(n_obs, mean = 0, sd = 1)          # diffusion shocks
  list(R = mu + sigma * eps + Y * Z, Y = Y, Z = Z)
}

#' The three parameter regimes used throughout (daily returns).
REGIMES <- list(
  low    = c(mu = 0.0003, sigma = 0.007, p = 0.01, mu_J = -0.015, sigma_J = 0.01),
  medium = c(mu = 0.0005, sigma = 0.010, p = 0.02, mu_J = -0.030, sigma_J = 0.02),
  high   = c(mu = 0.0000, sigma = 0.020, p = 0.06, mu_J = -0.060, sigma_J = 0.04)
)

PARAM_NAMES <- c("mu", "sigma", "p", "mu_J", "sigma_J")
