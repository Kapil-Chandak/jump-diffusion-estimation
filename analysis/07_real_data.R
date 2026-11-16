## Real data: JNJ, SPY and TSLA daily log-returns, 2015-2024 (as in the course report).
## 1. All four estimators on each stock.
## 2. Model check: are the days flagged as jumps spread evenly over the years, as independent jumps
##    would be? (The same test is validated on simulated data in analysis/06_independence_test.R.)
## Needs internet access (Yahoo Finance via quantmod); a few minutes because of Stan.
## Output: results/real_data.csv, results/independence_test_real.csv, results/flagged_days_per_year.csv,
##         figures/real_data_jump_prob.png
source("R/load.R")
library(quantmod)
tickers <- c("JNJ", "SPY", "TSLA")
px <- do.call(merge, lapply(tickers, function(s) Ad(getSymbols(s, from = "2015-01-01", to = "2024-12-31", auto.assign = FALSE))))
ret <- na.omit(diff(log(px))); colnames(ret) <- tickers
years <- format(index(ret), "%Y")

## 1. Estimates
rows <- list(); em_fit <- list()
for (s in tickers) {
  R <- as.numeric(ret[, s])
  em_fit[[s]] <- fit_jd_em(R)$theta
  set.seed(2025)
  st <- fit_jd_stan(R, chains = 4)
  est <- list(MLE = fit_jd_mle(R)$theta, EM = em_fit[[s]],
              Gibbs = fit_jd_gibbs(R)$summary[, "median"], Stan = st$summary[, "median"])
  for (m in names(est)) rows[[length(rows) + 1]] <- data.frame(ticker = s, method = m, t(est[[m]][PARAM_NAMES]),
                                                              stan_divergences = if (m == "Stan") st$divergences else NA)
}
out <- do.call(rbind, rows)
write.csv(out, "results/real_data.csv", row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)

## 2. Independence check
jump_prob <- lapply(tickers, function(s) jd_jump_prob(em_fit[[s]], as.numeric(ret[, s]))); names(jump_prob) <- tickers
per_year <- sapply(tickers, function(s) tapply(jump_prob[[s]] > 0.5, years, sum))
write.csv(data.frame(year = rownames(per_year), per_year), "results/flagged_days_per_year.csv", row.names = FALSE)
set.seed(707)
test <- do.call(rbind, lapply(tickers, function(s) {
  r <- jd_flag_test(jump_prob[[s]] > 0.5, years)
  data.frame(ticker = s, flagged = r[["flagged"]], p_value = r[["p_value"]])
}))
write.csv(test, "results/independence_test_real.csv", row.names = FALSE)
print(per_year); print(test, row.names = FALSE)

png("figures/real_data_jump_prob.png", width = 1800, height = 1400, res = 200)
par(mfrow = c(3, 1), mar = c(3, 4.5, 2.5, 1))
for (s in tickers) {
  plot(index(ret), jump_prob[[s]], type = "h", ylim = c(0, 1), xlab = "", ylab = "P(jump | return)",
       main = sprintf("%s: probability that each day contains a jump (EM estimates)", s))
  abline(h = 0.5, lty = 2, col = "grey40")
}
dev.off()
