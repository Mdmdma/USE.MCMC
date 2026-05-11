# Benchmark the MCMC sampling path across three implementation states:
#   baseline    -> current main, restart-style burn-in, pure R loop
#   rm-burnin   -> R loop with Robbins-Monro burn-in (no C++ yet)
#   rcpp-R      -> same as rm-burnin but on the rcpp branch (sanity check)
#   rcpp-cpp    -> C++ inner loop via Rcpp(Armadillo)
#
# Run from the package root:
#   Rscript development_scripts/benchmark_mcmc.R <git_label>
# e.g.
#   Rscript development_scripts/benchmark_mcmc.R baseline
#
# Appends rows to development_scripts/benchmark_results.csv.

suppressPackageStartupMessages({
  library(devtools)
  library(magrittr)
  library(sf)
  library(terra)
})

args <- commandArgs(trailingOnly = TRUE)
git_label <- if (length(args) >= 1) args[[1]] else "baseline"
n_iter <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
engine_for_label <- if (length(args) >= 3) args[[3]] else NA_character_

pkg_root <- getwd()
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  stop("Run this script from the package root (DESCRIPTION not found in ", pkg_root, ")")
}
devtools::load_all(pkg_root, quiet = TRUE)

git_sha <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE),
  error = function(e) "unknown"
)
results_csv <- file.path(pkg_root, "development_scripts", "benchmark_results.csv")

# ---- shared workload setup (runs once, outside any timed region) ----------

set.seed(42)

env.data.raster <- USE.MCMC::Worldclim_tmp %>%
  terra::rast(type = "xyz") %>%
  round(2)

rpc <- rastPCA(env.data.raster, stand = TRUE)
env.data.raster.with.pc <- c(env.data.raster, rpc$PCs)

env.with.pc.sf <- terra::as.data.frame(c(env.data.raster, rpc$PCs), xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y"))

set.seed(42)
virtual.presence.data <- getVirtualSpeciesPresencePoints(
  env.data = env.data.raster.with.pc, n.samples = 300
)
virtual.presence.points <- virtual.presence.data$sample.points

# Pre-fit GMMs + build density and proposal for the pure-chain benchmark.
# This mirrors the setup inside paSamplingMcmc() so the timed call only
# exercises the Metropolis inner loop.
fit_density_and_proposal <- function(dimensions) {
  env.with.pc.sf.subsampled <- env.with.pc.sf[
    stats::runif(min(nrow(env.with.pc.sf), 2000), 1, nrow(env.with.pc.sf)),
  ]
  env.data.cleaned <- sf::st_drop_geometry(env.with.pc.sf[dimensions])
  env.data.cleaned.subsampled <- sf::st_drop_geometry(
    env.with.pc.sf.subsampled[dimensions]
  )

  environmental.data.model <- mclust::densityMclust(
    env.data.cleaned.subsampled, plot = FALSE, verbose = FALSE
  )
  environmental.densities <- mclust::predict.densityMclust(
    environmental.data.model, env.data.cleaned
  )
  environment.threshold <- stats::quantile(environmental.densities, 0.001)

  virtual.presence.points.pc <- terra::extract(
    env.data.raster.with.pc, virtual.presence.points, bind = TRUE
  ) %>%
    sf::st_as_sf()
  species.model <- mclust::densityMclust(
    sf::st_drop_geometry(virtual.presence.points.pc[dimensions]),
    plot = FALSE, verbose = FALSE
  )
  species.cutoff.threshold <- stats::quantile(species.model$density, 0.95)

  densityFn <- mclustDensityFunction(
    env.model = environmental.data.model,
    species.model = species.model,
    dim = dimensions,
    threshold = environment.threshold,
    species.cutoff.threshold = species.cutoff.threshold
  )

  cov.mat <- 0.075 * stats::cov(
    sf::st_drop_geometry(env.with.pc.sf)[dimensions]
  )
  proposalFn <- addHighDimGaussian(
    cov.mat = cov.mat, dim = length(dimensions)
  )

  list(
    dataset = env.with.pc.sf,
    densityFunction = densityFn,
    proposalFunction = proposalFn,
    dimensions = dimensions
  )
}

# ---- timing helper ---------------------------------------------------------

time_repeated <- function(expr, n) {
  expr <- substitute(expr)
  caller <- parent.frame()
  times <- numeric(n)
  for (i in seq_len(n)) {
    times[i] <- system.time(eval(expr, caller), gcFirst = TRUE)[["elapsed"]]
  }
  list(median = stats::median(times), iqr = stats::IQR(times), n = n, all = times)
}

append_row <- function(metric, sub_case, engine, timing) {
  row <- data.frame(
    git_sha = git_sha,
    git_label = git_label,
    metric = metric,
    sub_case = sub_case,
    median_seconds = timing$median,
    iqr_seconds = timing$iqr,
    n_iter = timing$n,
    engine = engine,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  write_header <- !file.exists(results_csv)
  utils::write.table(
    row, results_csv,
    sep = ",", row.names = FALSE, col.names = write_header,
    append = !write_header, quote = TRUE
  )
  cat(sprintf(
    "[%s] %-22s %-14s median=%.3fs  iqr=%.3fs  (n=%d)\n",
    git_label, metric, sub_case, timing$median, timing$iqr, timing$n
  ))
}

# Detect whether mcmcSampling supports an `engine` argument (T2 onwards).
sampler_has_engine <- "engine" %in% names(formals(mcmcSampling))
default_engine <- if (sampler_has_engine) "auto" else NA_character_
engine_for_label <- if (!is.na(engine_for_label)) engine_for_label else default_engine

call_mcmc <- function(args) {
  if (sampler_has_engine && !is.na(engine_for_label)) {
    args$engine <- engine_for_label
  }
  do.call(mcmcSampling, args)
}

call_pa <- function(args) {
  if (sampler_has_engine && !is.na(engine_for_label) &&
      "engine" %in% names(formals(paSamplingMcmc))) {
    args$engine <- engine_for_label
  }
  do.call(paSamplingMcmc, args)
}

# ---- pure chain runtime ----------------------------------------------------

cat("=== Pure chain runtime ===\n")

for (dims in list(c("PC1", "PC2"), c("PC1", "PC2", "PC3"))) {
  setup <- fit_density_and_proposal(dims)
  args <- list(
    dataset = setup$dataset,
    dimensions = setup$dimensions,
    densityFunction = setup$densityFunction,
    proposalFunction = setup$proposalFunction,
    n.sample.points = 10000,
    burnIn = 1000,
    covariance.correction = 1,
    verbose = FALSE
  )
  # Warm up once (RNG batches, JIT, GMM precompute) then measure.
  invisible(call_mcmc(args))
  timing <- time_repeated(call_mcmc(args), n_iter)
  append_row(
    metric = "pure_chain",
    sub_case = sprintf("chain10k_d%d", length(dims)),
    engine = if (is.na(engine_for_label)) "R-legacy" else engine_for_label,
    timing = timing
  )
}

# ---- full workflow runtime -------------------------------------------------

cat("=== Full workflow runtime (paSamplingMcmc) ===\n")

full_args <- list(
  env.data.raster = env.data.raster,
  pres = virtual.presence.points,
  precomputed.pca = rpc,
  environmental.cutof.percentile = 0.001,
  num.chains = 4,
  num.cores = 1,
  chain.length = 10000,
  n.samples = 500,
  covariance.correction = 70,
  verbose = FALSE,
  dimensions = c("PC1", "PC2", "PC3")
)
invisible(call_pa(full_args))
timing <- time_repeated(call_pa(full_args), n_iter)
append_row(
  metric = "full_workflow",
  sub_case = "paSamplingMcmc_4chains_10k_d3",
  engine = if (is.na(engine_for_label)) "R-legacy" else engine_for_label,
  timing = timing
)

cat(sprintf("\nResults appended to %s\n", results_csv))

if (interactive()) {
  res <- utils::read.csv(results_csv, stringsAsFactors = FALSE)
  print(res)
}
