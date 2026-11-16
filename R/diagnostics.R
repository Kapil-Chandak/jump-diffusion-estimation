## Model check: are the jumps independent over time, as the model assumes?
## Days are flagged as jumps when P(Y_t = 1 | R_t) > 0.5 at the fitted parameters, and a chi-square test
## asks whether the flagged days are spread evenly across blocks of time (years). Under the model the
## flags are independent, so the counts per block are multinomial; with volatility clustering they bunch
## up in turbulent periods. The p-value is computed by simulation, which stays valid with few flags.

#' @param flags logical vector, one entry per day
#' @param block block label for each day (e.g. the calendar year)
jd_flag_test <- function(flags, block, B = 4999) {
  flags[is.na(flags)] <- FALSE
  counts <- tapply(flags, block, sum)
  n_days <- tapply(flags, block, length)
  if (sum(counts) < 2) return(c(flagged = sum(counts), p_value = NA))
  test <- suppressWarnings(chisq.test(counts, p = n_days / sum(n_days), simulate.p.value = TRUE, B = B))
  c(flagged = sum(counts), p_value = test$p.value)
}

#' Convenience wrapper: flag days at parameters theta and test them.
jd_independence_test <- function(theta, R, block, threshold = 0.5, B = 4999) {
  jd_flag_test(jd_jump_prob(theta, R) > threshold, block, B)
}
