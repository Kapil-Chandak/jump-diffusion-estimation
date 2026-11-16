// Jump-diffusion model, marginalised over the jump indicator (Stan has no discrete parameters).
// R_t ~ (1 - p) N(mu, sigma^2) + p N(mu + mu_J, sigma^2 + sigma_J^2)
// Priors are identical to the Gibbs sampler (R/gibbs.R), so the two Bayesian methods are directly comparable.
data {
  int<lower=1> T;
  vector[T] R;
}
parameters {
  real mu;
  real<lower=0> sigma2;
  real<lower=0, upper=1> p;
  real mu_J;
  real<lower=0> sigma_J2;
}
model {
  vector[T] lp0;
  vector[T] lp1;
  vector[T] m;
  mu ~ normal(0, 0.05);
  sigma2 ~ inv_gamma(2, 2e-4);
  p ~ beta(1, 9);
  mu_J ~ normal(0, 0.05);
  sigma_J2 ~ inv_gamma(2, 2e-4);
  lp0 = log1m(p) - 0.5 * log(2 * pi() * sigma2) - square(R - mu) / (2 * sigma2);
  lp1 = log(p) - 0.5 * log(2 * pi() * (sigma2 + sigma_J2)) - square(R - mu - mu_J) / (2 * (sigma2 + sigma_J2));
  m = fmax(lp0, lp1);                                // per-day log-sum-exp, vectorised
  target += sum(m + log(exp(lp0 - m) + exp(lp1 - m)));
}
generated quantities {
  real sigma = sqrt(sigma2);
  real sigma_J = sqrt(sigma_J2);
}
