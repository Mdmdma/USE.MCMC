#' MCMC sampling from a given dataset
#'
#' @param dataset sf dataframe from which the points are sampled
#' @param dimensions string vector containing the dimensions that should be included in the random walk
#' @param densityFunction Function that can take a point given as a numeric vector as input and returns the target density at that location.
#' @param proposalFunction Function that can take a point given as a numeric vector and a covariance adjuster as input and returns a new proposed point as a numeric vector.
#' @param n.sample.points Number of points to be sampled
#' @param burnIn Integer, sets the number of samples per adaptive burn in step. If set to 0, burn in is skipped
#' @param verbose Boolean to toggle progress updates
#' @param covariance.correction Integer, initial value of the covariance correction.
#' @param max.burnin.cycles Integer, maximum number of burn-in adaptation cycles before stopping with a warning. Prevents infinite loops when the target acceptance rate cannot be reached.
#' @returns A data.frame containing the sampled points with dimension columns and a density column
#' @export
#'
mcmcSampling <- function(dataset = NULL, dimensions= list(""), densityFunction = alwaysOne,
                         proposalFunction = addHighDimGaussian(dim = length(dimensions)),
                         n.sample.points = 0, burnIn = 1000, verbose = TRUE,
                         covariance.correction = 1, max.burnin.cycles = 50){
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
  if (!is.numeric(max.burnin.cycles) || length(max.burnin.cycles) != 1 || max.burnin.cycles < 1) {
    stop(paste0("'max.burnin.cycles' must be a positive integer, got ", deparse(max.burnin.cycles)), call. = FALSE)
  }

  n.dim <- length(dimensions)
  starting.index <- sample(nrow(dataset), 1)
  start.row <- sf::st_drop_geometry(dataset[starting.index, ])
  current.point <- as.numeric(start.row[, dimensions])
  names(current.point) <- dimensions
  current.density <- densityFunction(current.point)

  # burn in
  if(burnIn > 0) {
    cat("Burn in\n")
    points.accepted <- 0
    burnin.cycle <- 0
    if (verbose){
      pb.burnin <- utils::txtProgressBar(min = 0, max = burnIn, style = 3)
    }
   while (points.accepted / burnIn < 0.21 | points.accepted / burnIn > 0.25) {
      # the numbers of the condition depend on the exact threshold, Gelman, Roberts, and Gilks (1996) proposes 0.23 was optimal
      burnin.cycle <- burnin.cycle + 1
      if (burnin.cycle > max.burnin.cycles) {
        warning(paste0("Burn-in did not converge after ", max.burnin.cycles,
                       " cycles (last acceptance rate: ", round(points.accepted / burnIn, 3),
                       "). Proceeding with current covariance.correction = ",
                       round(covariance.correction, 6)), call. = FALSE)
        break
      }

      if (verbose){
            cat("\rThe current acceptance rate is", points.accepted /burnIn, ", the currend covariance adjustment factor is ", covariance.correction, "\n")
      }
      points.accepted <- 0
      for (i in 1:burnIn) {
        proposed.point <- proposalFunction(current.point, covariance.adjuster = covariance.correction, dim = dimensions)
        proposed.density <- densityFunction(proposed.point)
        if (acceptNextPoint(current.density, proposed.density)){
          current.point <- proposed.point
          current.density <- proposed.density
          points.accepted <- points.accepted + 1
        }

        if (verbose){
          cat("\rCurrent acceptance ratio: ", points.accepted / i)
          utils::setTxtProgressBar(pb.burnin, i)
        }
      }
      if (points.accepted / burnIn < 0.21) covariance.correction <- max(1e-10, covariance.correction * stats::rnorm(1, mean = 0.7, sd = 0.1))
      if (points.accepted / burnIn > 0.25) covariance.correction <- covariance.correction * stats::rnorm(1, mean = 1.3, sd = 0.1)
    }
  }
  else cat("Burn in skipped")

  cat("\nThe final covariance correction is ", covariance.correction, "\n")
  # Pre-allocate output matrix (dimensions + density)
  n.cols <- n.dim + 1
  sampled.points <- matrix(NA_real_, nrow = n.sample.points, ncol = n.cols)
  points.rejected <- 0
  iteration <- 0

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n.sample.points, style = 3)
  }

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
    if (verbose){
      cat("\rPoints rejected:", points.rejected, "Points sampled:", iteration)
      utils::setTxtProgressBar(pb, iteration)
    }
  }
  cat("\nPoints rejected: ", points.rejected)

  # Convert output matrix to data.frame
  sampled.points <- as.data.frame(sampled.points)
  colnames(sampled.points) <- c(dimensions, "density")
  return(sampled.points)
}
