#' Helper to create a Density function that uses mclust Gaussian mixtures
#'
#' @param env.model mclust gaussian mixture that uses points
#' @param species.model mclust gaussian mixture that uses points
#' @param dim string vector specifing the names of the dimensions
#' @param threshold sets the cutoff density from the environment
#' @param species.cutoff.threshold set the scaling factor with which the species model gets scaled before subtraction.
#'
#' @returns Function that can calculates the density at a point
#' @export
#'
mclustDensityFunction <- function(env.model = NULL, species.model = NULL, dim = "", threshold = 0.01, species.cutoff.threshold = 0.1){
  # Input validation
  if (is.null(env.model)) {
    stop("'env.model' must be provided (got NULL)", call. = FALSE)
  }
  if (!inherits(env.model, "densityMclust")) {
    stop(paste0("'env.model' must be a densityMclust object (from mclust::densityMclust), got '",
                paste(class(env.model), collapse = "/"), "'"), call. = FALSE)
  }
  if (is.null(species.model)) {
    stop("'species.model' must be provided (got NULL)", call. = FALSE)
  }
  if (!inherits(species.model, "densityMclust")) {
    stop(paste0("'species.model' must be a densityMclust object (from mclust::densityMclust), got '",
                paste(class(species.model), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.character(dim) || length(dim) < 1 || all(dim == "")) {
    stop("'dim' must be a character vector of dimension names", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0) {
    stop(paste0("'threshold' must be a positive number, got ", deparse(threshold)), call. = FALSE)
  }
  if (!is.numeric(species.cutoff.threshold) || length(species.cutoff.threshold) != 1 || species.cutoff.threshold <= 0) {
    stop(paste0("'species.cutoff.threshold' must be a positive number, got ", deparse(species.cutoff.threshold)), call. = FALSE)
  }

  # Precompute GMM parameters for fast evaluation
  env.pre <- precompute_gmm_params(env.model)
  species.pre <- precompute_gmm_params(species.model)
  threshold.floor <- threshold / 1000

  densityAtPointEstimator <- function(point){
    density <- fast_gmm_density(point, env.pre)
    if (density < threshold) return(threshold.floor)
    return(max(threshold.floor, 1 - fast_gmm_density(point, species.pre) / species.cutoff.threshold))
    }
  return(densityAtPointEstimator)
}

#' Precompute GMM parameters for fast density evaluation
#'
#' Extracts and precomputes inverse covariance matrices and log normalization
#' constants from a densityMclust model to avoid recomputing them on every call.
#'
#' @param model A densityMclust object
#' @returns A list with precomputed parameters
#' @keywords internal
precompute_gmm_params <- function(model) {
  G <- model$G
  d <- model$d
  pro <- model$parameters$pro
  mu <- model$parameters$mean
  sigma <- model$parameters$variance$sigma

  # Ensure sigma is a 3D array even for G=1
  if (length(base::dim(sigma)) == 2) {
    sigma <- array(sigma, dim = c(d, d, 1))
  }
  # Ensure mu is a matrix even for G=1
  if (is.null(base::dim(mu))) {
    mu <- matrix(mu, ncol = 1)
  }

  inv_sigma <- array(0, dim = c(d, d, G))
  log_norm <- numeric(G)

  for (k in seq_len(G)) {
    sig_k <- sigma[, , k]
    inv_sigma[, , k] <- solve(sig_k)
    log_norm[k] <- log(pro[k]) - 0.5 * d * log(2 * pi) - 0.5 * log(det(sig_k))
  }

  list(inv_sigma = inv_sigma, log_norm = log_norm, mean = mu, G = G, d = d)
}

#' Fast GMM density evaluation using precomputed parameters
#'
#' Computes the density of a Gaussian mixture model at a point using
#' precomputed inverse covariance matrices and log normalization constants.
#' Uses log-sum-exp for numerical stability.
#'
#' @param x Numeric vector, the point at which to evaluate the density.
#' @param pre List of precomputed parameters from \code{precompute_gmm_params}.
#' @returns Numeric scalar, the density at x.
#' @keywords internal
fast_gmm_density <- function(x, pre) {
  log_densities <- numeric(pre$G)
  for (k in seq_len(pre$G)) {
    diff <- x - pre$mean[, k]
    log_densities[k] <- pre$log_norm[k] - 0.5 * sum(diff * (pre$inv_sigma[, , k] %*% diff))
  }
  max_log <- max(log_densities)
  exp(max_log + log(sum(exp(log_densities - max_log))))
}

#' Batched GMM density evaluation using precomputed parameters
#'
#' Vectorized variant of \code{fast_gmm_density} for evaluating many points at
#' once. Operates per component on the whole batch of points, then combines via
#' log-sum-exp. Substantially faster than looping single-point evaluations when
#' the caller already has all query points in hand (e.g. plotting density rasters,
#' threshold sweeps, or sweep-style diagnostics in the vignettes).
#'
#' @param X Numeric matrix (n x d) or data.frame coercible to one — query points
#'   in rows, dimensions in columns.
#' @param pre List of precomputed parameters from \code{precompute_gmm_params}.
#' @returns Numeric vector of length n with the density at each row of X.
#' @keywords internal
fast_gmm_density_batch <- function(X, pre) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) X <- matrix(X, ncol = pre$d)
  n <- nrow(X)
  log_densities <- matrix(NA_real_, nrow = n, ncol = pre$G)
  for (k in seq_len(pre$G)) {
    diffs <- sweep(X, 2L, pre$mean[, k], FUN = "-")              # n x d
    # quadratic form: rowSums(diffs * (diffs %*% inv_sigma)) gives diff^T A diff per row.
    log_densities[, k] <- pre$log_norm[k] - 0.5 * rowSums(diffs * (diffs %*% pre$inv_sigma[, , k]))
  }
  max_log <- do.call(pmax, as.data.frame(log_densities))
  exp(max_log + log(rowSums(exp(log_densities - max_log))))
}
