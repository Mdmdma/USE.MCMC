#' USE.MCMC: pseudo-absence sampling in environmental space
#'
#' Forks the original USE package to support pseudo-absence sampling in
#' higher-dimensional environmental spaces. Two backends are exposed: a
#' Metropolis MCMC sampler (`paSamplingMcmc()`, with the underlying
#' `mcmcSampling()` engine) and a nearest-neighbor sampler (`paSamplingNn()`).
#'
#' The MCMC inner loop is implemented in C++ via Rcpp + RcppArmadillo and is
#' dispatched automatically whenever the density and proposal closures are
#' built by the package's own factories (`mclustDensityFunction()` and
#' `addHighDimGaussian()`). User-supplied closures continue to run on the R
#' reference loop. Engine selection is controlled by the `engine` argument
#' and defaults to `"auto"`.
#'
#' @useDynLib USE.MCMC, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
