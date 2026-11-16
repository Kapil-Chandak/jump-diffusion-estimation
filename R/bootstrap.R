## Nonparametric bootstrap standard errors for a point estimator.
## Days are resampled with replacement (valid here because the model treats days as independent),
## the estimator is refitted to each resample, and the standard deviation of the estimates is the SE.

#' @param R       daily returns
#' @param fit_fun estimator returning a list with element `theta` (e.g. fit_jd_mle, fit_jd_em)
#' @param B       number of bootstrap resamples
jd_bootstrap <- function(R, fit_fun, B = 300) {
  est <- t(replicate(B, {
    th <- tryCatch(fit_fun(sample(R, replace = TRUE))$theta, error = function(e) rep(NA_real_, 5))
    unname(th[PARAM_NAMES])
  }))
  colnames(est) <- PARAM_NAMES
  list(se = apply(est, 2, sd, na.rm = TRUE), estimates = est, failures = sum(!complete.cases(est)))
}
