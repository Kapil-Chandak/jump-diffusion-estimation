## Load all functions of the project. Usage: source("R/load.R") from the repository root.
.repo_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "load.R")) && dirname(d) != d) d <- dirname(d)
  d
}
for (f in c("simulate.R", "likelihood.R", "mle.R", "em.R", "gibbs.R", "stan.R", "bootstrap.R", "diagnostics.R"))
  source(file.path(.repo_root(), "R", f))

#' Parallel lapply: forked workers on macOS/Linux; on Windows (no forking) it runs sequentially.
par_map <- function(X, FUN, cores = 1) {
  if (cores > 1 && .Platform$OS.type != "windows")
    parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
  else lapply(X, FUN)
}
