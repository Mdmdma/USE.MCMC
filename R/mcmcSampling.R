# Minimum proposal-scale floor; keep in sync with MIN_COV_CORRECTION in src/mcmc_loop.cpp.
MIN_COV_CORRECTION <- 1e-10

#' MCMC sampling from a given dataset
#'
#' @param dataset sf dataframe from which the points are sampled
#' @param dimensions string vector containing the dimensions that should be included in the random walk
#' @param densityFunction Function that can take a point given as a numeric vector as input and returns the target density at that location.
#' @param proposalFunction Function that can take a point given as a numeric vector and a covariance adjuster as input and returns a new proposed point as a numeric vector.
#' @param n.sample.points Number of points to be sampled
#' @param burnIn Integer, number of Robbins-Monro burn-in adaptation steps performed before sampling. During each step the proposal scale is adjusted toward target acceptance 0.234 (Roberts/Rosenthal 2009). Set to 0 to skip adaptation and start sampling immediately at the user-supplied `covariance.correction`.
#' @param verbose Boolean to toggle progress updates
#' @param covariance.correction Integer, initial value of the covariance correction.
#' @param max.burnin.cycles Deprecated. Retained for backwards compatibility; ignored by the current Robbins-Monro burn-in.
#' @param engine One of `"auto"` (default), `"R"`, or `"cpp"`. `"auto"` picks the C++ inner loop when both `densityFunction` and `proposalFunction` are built by `mclustDensityFunction()` and `addHighDimGaussian()` (they carry the required `rcpp_spec` attribute) and falls back to the R loop otherwise. `"cpp"` forces the C++ path and errors if a custom closure is supplied. `"R"` forces the pure-R reference loop.
#' @returns A data.frame containing the sampled points with dimension columns and a density column
#' @export
mcmcSampling <- function(dataset = NULL, dimensions= list(""), densityFunction = alwaysOne,
                         proposalFunction = addHighDimGaussian(dim = length(dimensions)),
                         n.sample.points = 0, burnIn = 1000, verbose = TRUE,
                         covariance.correction = 1, max.burnin.cycles = 50,
                         engine = c("auto", "R", "cpp")){
  # Input validation
  if (is.null(dataset)) {
    stop("'dataset' must be provided (got NULL)", call. = FALSE)
  }
  if (!is.data.frame(dataset) && !inherits(dataset, "sf")) {
    stop(paste0("'dataset' must be a data.frame or sf object, got '",
                paste(class(dataset), collapse = "/"), "'"), call. = FALSE)
  }
  if (nrow(dataset) == 0) {
    stop("'dataset' must have at least one row", call. = FALSE)
  }
  if (!is.character(dimensions) && !is.list(dimensions)) {
    stop(paste0("'dimensions' must be a character vector, got '",
                paste(class(dimensions), collapse = "/"), "'"), call. = FALSE)
  }
  dimensions <- unlist(dimensions)
  if (!is.character(dimensions) || length(dimensions) < 1 || all(dimensions == "")) {
    stop("'dimensions' must be a non-empty character vector of column names", call. = FALSE)
  }
  check_columns_exist(dataset, dimensions, "dataset", "dimensions")
  if (!is.function(densityFunction)) {
    stop(paste0("'densityFunction' must be a function, got '",
                paste(class(densityFunction), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.function(proposalFunction)) {
    stop(paste0("'proposalFunction' must be a function, got '",
                paste(class(proposalFunction), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.numeric(n.sample.points) || length(n.sample.points) != 1 || n.sample.points < 1) {
    stop(paste0("'n.sample.points' must be a positive number, got ", deparse(n.sample.points)), call. = FALSE)
  }
  if (is.logical(burnIn)) {
    stop("'burnIn' must be an integer, not logical (TRUE is coerced to 1 which causes issues)", call. = FALSE)
  }
  if (!is.numeric(burnIn) || length(burnIn) != 1 || burnIn < 0) {
    stop(paste0("'burnIn' must be a non-negative number, got ", deparse(burnIn)), call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop(paste0("'verbose' must be a single logical value, got '",
                paste(class(verbose), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.numeric(covariance.correction) || length(covariance.correction) != 1 || covariance.correction <= 0) {
    stop(paste0("'covariance.correction' must be a positive number, got ", deparse(covariance.correction)), call. = FALSE)
  }
  if (!missing(max.burnin.cycles) && !identical(max.burnin.cycles, 50)) {
    warning("'max.burnin.cycles' is deprecated and ignored; burn-in now uses single-pass Robbins-Monro adaptation.", call. = FALSE)
  }
  engine <- match.arg(engine)

  density.spec <- attr(densityFunction, "rcpp_spec")
  proposal.spec <- attr(proposalFunction, "rcpp_spec")
  cpp.eligible <- !is.null(density.spec) &&
                  !is.null(proposal.spec) &&
                  identical(density.spec$type, "mclust_density") &&
                  identical(proposal.spec$type, "gaussian_proposal")

  if (engine == "cpp" && !cpp.eligible) {
    stop("engine = 'cpp' requires densityFunction from mclustDensityFunction() and proposalFunction from addHighDimGaussian(). Use engine = 'R' or 'auto' for custom closures.", call. = FALSE)
  }
  if (engine == "auto" && cpp.eligible) {
    engine <- "cpp"
  } else if (engine == "auto") {
    engine <- "R"
  }

  n.dim <- length(dimensions)
  starting.index <- sample(nrow(dataset), 1)
  start.row <- sf::st_drop_geometry(dataset[starting.index, ])
  current.point <- as.numeric(start.row[, dimensions])
  names(current.point) <- dimensions

  if (engine == "cpp") {
    if (verbose) cat("Running C++ inner loop (Rcpp)\n")
    tryCatch(
      chol(proposal.spec$cov),
      error = function(e) {
        stop("proposal covariance must be positive-definite for engine = 'cpp'; ",
             "got a matrix that failed Cholesky decomposition. ",
             "Check the 'cov.mat' argument to addHighDimGaussian().",
             call. = FALSE)
      }
    )
    res <- mcmc_loop_cpp(
      start_point = current.point,
      n_sample_points = as.integer(n.sample.points),
      burn_in = as.integer(burnIn),
      covariance_correction = covariance.correction,
      env_inv_sigma = density.spec$env$inv_sigma,
      env_log_norm = density.spec$env$log_norm,
      env_means = density.spec$env$mean,
      env_threshold = density.spec$threshold,
      sp_inv_sigma = density.spec$species$inv_sigma,
      sp_log_norm = density.spec$species$log_norm,
      sp_means = density.spec$species$mean,
      sp_cutoff = density.spec$species_cutoff,
      floor_value = density.spec$floor,
      proposal_mean = as.numeric(proposal.spec$mean),
      proposal_cov = proposal.spec$cov
    )
    sampled.points <- as.data.frame(res$samples)
    colnames(sampled.points) <- c(dimensions, "density")
    if (verbose) {
      cat(sprintf("\nFinal covariance.correction: %g; sampling rejected: %d\n",
                  res$covariance_correction, res$sampling_rejected))
    }
    return(sampled.points)
  }

  current.density <- densityFunction(current.point)

  # Burn-in: single-pass Robbins-Monro adaptation on log(covariance.correction).
  # After each step, log(c) += gamma_t * (accept - target), gamma_t = 1/(t+1)^0.6.
  # Target acceptance 0.234 (Roberts/Rosenthal 2009 optimal for high-d random-walk MH).
  if(burnIn > 0) {
    if (verbose) cat("Burn in (Robbins-Monro)\n")
    target.acceptance <- 0.234
    rm.exponent <- 0.6
    log.cov.correction <- log(covariance.correction)
    points.accepted <- 0
    if (verbose){
      pb.burnin <- utils::txtProgressBar(min = 0, max = burnIn, style = 3)
    }
    progress.stride <- max(1L, as.integer(burnIn %/% 100L))
    for (t in seq_len(burnIn)) {
      proposed.point <- proposalFunction(current.point, covariance.adjuster = exp(log.cov.correction), dim = dimensions)
      proposed.density <- densityFunction(proposed.point)
      accepted <- acceptNextPoint(current.density, proposed.density)
      if (accepted){
        current.point <- proposed.point
        current.density <- proposed.density
        points.accepted <- points.accepted + 1
      }
      gamma_t <- 1 / (t + 1)^rm.exponent
      log.cov.correction <- log.cov.correction + gamma_t * ((if (accepted) 1 else 0) - target.acceptance)

      if (verbose && (t %% progress.stride == 0L || t == burnIn)){
        cat("\rCurrent acceptance ratio: ", points.accepted / t, " covariance.correction: ", exp(log.cov.correction))
        utils::setTxtProgressBar(pb.burnin, t)
      }
    }
    covariance.correction <- max(MIN_COV_CORRECTION, exp(log.cov.correction))
  }
  else if (verbose) cat("Burn in skipped")

  if (verbose) cat("\nThe final covariance correction is ", covariance.correction, "\n")
  # Pre-allocate output matrix (dimensions + density)
  n.cols <- n.dim + 1
  sampled.points <- matrix(NA_real_, nrow = n.sample.points, ncol = n.cols)
  points.rejected <- 0
  iteration <- 0

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n.sample.points, style = 3)
  }
  sampling.stride <- max(1L, as.integer(n.sample.points %/% 100L))

  while (iteration < n.sample.points) {
    proposed.point <- proposalFunction(current.point, covariance.adjuster = covariance.correction, dim = dimensions)
    proposed.density <- densityFunction(proposed.point)
    iteration <- iteration + 1
    if (acceptNextPoint(current.density, proposed.density)){
      current.point <- proposed.point
      current.density <- proposed.density
    }
    else {
      points.rejected <- points.rejected + 1
    }
    sampled.points[iteration, ] <- c(current.point, current.density)
    if (verbose && (iteration %% sampling.stride == 0L || iteration == n.sample.points)){
      cat("\rPoints rejected:", points.rejected, "Points sampled:", iteration)
      utils::setTxtProgressBar(pb, iteration)
    }
  }
  if (verbose) cat("\nPoints rejected: ", points.rejected)

  # Convert output matrix to data.frame
  sampled.points <- as.data.frame(sampled.points)
  colnames(sampled.points) <- c(dimensions, "density")
  return(sampled.points)
}
