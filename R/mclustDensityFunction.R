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

  densityAtPointEstimator <- function(point){
    point <- as.matrix(sf::st_drop_geometry(point[dim]))
    density <- mclust::predict.densityMclust(env.model, point)
    if (density < threshold ) return(threshold / 1000)
    return(max(threshold / 1000, 1 - mclust::predict.densityMclust(species.model, point)  / species.cutoff.threshold)) # this approach is dubeaous
    }
  return(densityAtPointEstimator)
}
