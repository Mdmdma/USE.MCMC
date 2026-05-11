# USE.MCMC 0.0.3

## MCMC backend
* New C++ (Rcpp + RcppArmadillo) inner loop with auto-dispatch via a new
  `engine = c("auto", "R", "cpp")` argument on `mcmcSampling()` and
  `paSamplingMcmc()`. The default `"auto"` picks the C++ path whenever both
  density and proposal closures are built by `mclustDensityFunction()` and
  `addHighDimGaussian()`; custom user closures continue to use the R loop.
* Burn-in uses single-pass Robbins-Monro adaptation (target acceptance 0.234).
* `paSamplingMcmc` setup builds the env+PC frame via `terra::as.data.frame`
  on a multi-layer SpatRaster, dropping the previous `sf::st_join` hop.
* Verbose output in `mcmcSampling()` is gated by `verbose` and throttled to
  ~100 updates so parallel runs no longer serialize tens of thousands of
  `cat()` calls.
* Batched `fast_gmm_density_batch()` helper in `mclustDensityFunction.R` for
  vectorized GMM density evaluation.

## NN backend
* Fixes the previously crashing default branch
  (`nn.based.presence.exclusion = TRUE`) of `paSamplingNn()`.
* Fixes `step.x`/`step.y` collapsing y-noise, `distance.threshold` being
  unconditionally overwritten, hardcoded `c("PC1", "PC2")` ignoring the
  `dimensions` parameter, and an off-by-one in `maxResNn()`'s
  neighbor-distance subset.
* `optimalDistanceThresholdNn()` uses `index.for.cutof` consistently and a
  partial sort instead of a full sort.
* `paSamplingNn()` drops NA presence rows before `FNN::get.knnx()` and errors
  cleanly when fewer than two non-NA points remain.

## Cleanups and tooling
* Drops four unused R files
  (`gaussianMixureDensityFunction.R`, `analyseAreaInHighDim.R`,
  `mapBackOnRealPoints.R`, `EuclidianMetric.R`) and seven superseded
  development scripts.
* Cluster precomputation scripts modernized: integer `burnIn`,
  `terra::as.data.frame` setup, `USE_MCMC_CHAINS_DIR` env var, CLI options
  matching current library defaults.
* Methods-overview compound figure added to the MCMC vignette;
  vignette gifs ship as `knitr::include_graphics()` resources so
  `devtools::build_vignettes()` copies them next to the rendered HTML.

# USE.MCMC 0.0.1
* First release.
