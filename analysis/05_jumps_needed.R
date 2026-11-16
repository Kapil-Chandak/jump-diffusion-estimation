## How much data does it take? Accuracy of the jump frequency p against the expected number of
## jumps in the sample (T * p), for small and large jumps relative to ordinary volatility.
## EM only (it reaches the same maximum as the MLE and is fast).
## Usage: Rscript analysis/05_jumps_needed.R [reps=100] [cores=2]. Output: results/jumps_needed.csv,
## figures/jumps_needed.png
source("R/load.R")
grid <- expand.grid(T = c(500, 1000, 2000, 5000, 10000), p = c(0.005, 0.01, 0.02, 0.05), size = c(2, 4))
grid$sigma <- 0.01; grid$mu_J <- -grid$size * grid$sigma; grid$sigma_J <- grid$sigma
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) > 0) as.integer(args[1]) else 100
CORES <- if (length(args) > 1) as.integer(args[2]) else 2
cell <- function(i) {
  g <- grid[i, ]
  set.seed(30300 + i)                          # one seed per cell: same results on any number of cores
  err <- replicate(REPS, {
    R <- simulate_jump_diffusion(g$T, 0, g$sigma, g$p, g$mu_J, g$sigma_J)$R
    fit_jd_em(R)$theta[["p"]] / g$p - 1 })
  data.frame(g, expected_jumps = g$T * g$p, median_abs_rel_error = median(abs(err)),
             share_within_50pct = mean(abs(err) < 0.5))
}
res <- do.call(rbind, par_map(seq_len(nrow(grid)), cell, CORES))
write.csv(res, "results/jumps_needed.csv", row.names = FALSE)
png("figures/jumps_needed.png", width = 1600, height = 1000, res = 220)
par(mar = c(4.5, 4.5, 2.5, 1))
plot(NA, xlim = range(res$expected_jumps), ylim = c(0, 1), log = "x", xlab = "Expected number of jumps in the sample (T x p)",
     ylab = "Share of estimates of p within 50% of the truth", main = "Jump frequency needs enough jumps, and large enough jumps")
cols <- c(`2` = "#c0392b", `4` = "#2471a3")
for (s in c(2, 4)) { d <- res[res$size == s, ]; d <- d[order(d$expected_jumps), ]
  points(d$expected_jumps, d$share_within_50pct, pch = 19, col = cols[as.character(s)])
  lo <- lowess(log(d$expected_jumps), d$share_within_50pct, f = 0.6); lines(exp(lo$x), pmin(pmax(lo$y, 0), 1), col = cols[as.character(s)], lwd = 2) }
legend("bottomright", c("jumps 2x ordinary daily volatility", "jumps 4x ordinary daily volatility"), col = cols, pch = 19, lwd = 2, bty = "n")
dev.off()
print(res[order(res$size, res$expected_jumps), c("size", "T", "p", "expected_jumps", "median_abs_rel_error", "share_within_50pct")], digits = 3, row.names = FALSE)
