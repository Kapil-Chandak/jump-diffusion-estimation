## Run everything, in order, from the repository root:   Rscript run_all.R
## Approximate times on a laptop with 4 cores are given for each step (on Windows the parallel steps run on
## one core and take longer). Outputs go to results/ and figures/; a log is written to results/run_all.log.
## Work from the repository root wherever this file is run from (Rscript, source() or RStudio)
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(f)) f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!length(f) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    f <- rstudioapi::getSourceEditorContext()$path
  if (length(f) && nzchar(f)) setwd(dirname(normalizePath(f)))
  if (!file.exists("R/load.R")) stop("Run this file from the jump-diffusion-estimation folder (could not find R/load.R).")
})
CORES <- max(1, parallel::detectCores() - 1)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
log_file <- file("results/run_all.log", open = "wt"); sink(log_file, split = TRUE); sink(log_file, type = "message")
step <- function(name, expr) {
  cat("\n==========", name, "==========\n"); t0 <- Sys.time()
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("---------- %s: %s in %.1f min\n", name, if (ok) "done" else "FAILED", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
run <- function(script, args = character()) {
  status <- system2("Rscript", c(script, args), stdout = "", stderr = "")
  if (status != 0) stop(script, " exited with status ", status)
}
step("0. tests (about 3 min)",                          run("tests/run_tests.R"))
step("1. single dataset per regime (about 5 min)",      run("analysis/01_single_dataset.R"))
step("2. bootstrap standard errors (about 15 min)",     run("analysis/02_bootstrap.R", "300"))
step("3. prior sensitivity (about 3 min)",              run("analysis/03_prior_sensitivity.R"))
step("4. Monte Carlo study (about 30-60 min)",          run("analysis/04_monte_carlo.R", c("200", "25", CORES)))
step("5. jumps needed (about 10 min)",                  run("analysis/05_jumps_needed.R", c("100", CORES)))
step("6. independence test, simulated (about 5 min)",   run("analysis/06_independence_test.R", c("200", CORES)))
step("7. real data (needs internet; about 10 min)",     run("analysis/07_real_data.R"))
step("8. summary table and figure (seconds)",           run("analysis/08_summarise.R"))
sink(type = "message"); sink()
cat("\nAll steps attempted. See results/run_all.log for any failures.\n")
