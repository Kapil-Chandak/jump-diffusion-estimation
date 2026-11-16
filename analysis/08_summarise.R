## Summary table and figure from the Monte Carlo study (results/monte_carlo.rds).
## Output: results/monte_carlo_summary.csv, figures/monte_carlo_p.png
source("R/load.R")
args <- commandArgs(trailingOnly = TRUE)
x <- readRDS(if (length(args)) args[1] else "results/monte_carlo.rds")
methods <- c("report_MLE", "report_EM", "MLE", "EM", "Gibbs", "Stan")
x <- x[x$method %in% methods, ]
x$ratio <- x$estimate / x$truth
x$covered <- ifelse(is.na(x$lo), NA, x$lo <= x$truth & x$truth <= x$hi)

tab <- do.call(rbind, lapply(names(REGIMES), function(reg) do.call(rbind, lapply(methods, function(m) {
  d <- x[x$regime == reg & x$method == m, ]
  if (!nrow(d)) return(NULL)
  g <- function(par) d[d$parameter == par, ]
  data.frame(regime = reg, method = m, datasets = nrow(g("p")),
             median_p_ratio          = median(g("p")$ratio),
             median_abs_error_p      = median(abs(g("p")$ratio - 1)),
             median_muJ_ratio        = median(g("mu_J")$ratio),
             median_sigmaJ_ratio     = median(g("sigma_J")$ratio),
             share_sigmaJ_collapsed  = mean(g("sigma_J")$ratio < 0.1),
             coverage_95_p           = if (all(is.na(g("p")$covered))) NA else mean(g("p")$covered),
             coverage_95_sigmaJ      = if (all(is.na(g("sigma_J")$covered))) NA else mean(g("sigma_J")$covered))
}))))
write.csv(tab, "results/monte_carlo_summary.csv", row.names = FALSE)
print(format(tab, digits = 3), row.names = FALSE)

## MLE and EM maximise the same likelihood: how often do they reach the same maximum?
ll <- merge(x[x$method == "MLE" & x$parameter == "p", c("rep", "regime", "loglik")],
            x[x$method == "EM" & x$parameter == "p", c("rep", "regime", "loglik")], by = c("rep", "regime"))
cat(sprintf("MLE and EM log-likelihoods within 0.01: %.1f%% of datasets\n", 100 * mean(abs(ll$loglik.x - ll$loglik.y) < 0.01)))
cat("Stan divergences (total over fits):", sum(x$divergences[x$method == "Stan" & x$parameter == "p"], na.rm = TRUE),
    "in", sum(x$method == "Stan" & x$parameter == "p"), "fits\n")

## Figure: estimated / true jump frequency p, by method and regime (log scale)
x$method <- factor(x$method, methods)
png("figures/monte_carlo_p.png", width = 1800, height = 900, res = 200)
par(mfrow = c(1, 3), mar = c(7.5, 4.5, 3, 1))
lab <- c(report_MLE = "Original MLE", report_EM = "Original EM", MLE = "MLE", EM = "EM", Gibbs = "Gibbs", Stan = "Stan")
col <- c("#c0392b", "#e59866", "#2471a3", "#5499c7", "#1e8449", "#52be80")
lims <- c(0.02, 50)                     # shared axis; estimates outside it are drawn at the edge
for (r in names(REGIMES)) {
  d <- droplevels(x[x$regime == r & x$parameter == "p", ])
  d$ratio <- pmin(pmax(d$ratio, lims[1]), lims[2])
  boxplot(ratio ~ method, d, log = "y", las = 2, ylim = lims, yaxt = "n", names = lab[levels(d$method)], col = col[methods %in% levels(d$method)],
          outline = FALSE, ylab = "estimated p / true p", main = sprintf("%s regime (true p = %s)", r, REGIMES[[r]][["p"]]), xlab = "")
  axis(2, at = c(0.02, 0.1, 0.5, 1, 2, 10, 50), labels = c("<0.02", "0.1", "0.5", "1", "2", "10", ">50"), las = 1)
  abline(h = 1, lty = 2)
}
dev.off()
