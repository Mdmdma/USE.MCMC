#' addHighDimGaussian
#'
#' Adds high dimensional Gaussian noise to specific parameters of a point given as a numeric vector.
#'
#'
#' @param dim integer that specifies the number of dimensions that should be altered
#' @param mean.vec vector of the means of the Gaussian to be added
#' @param cov.mat covariance matrix of the Gaussian to be added
#' @param batch.size integer, number of random vectors to pre-generate at a time for efficiency
#'
#' @returns function that takes a point given as a numeric vector and returns it with Gaussian noise added
#' @export
#'
addHighDimGaussian <- function(dim = 0, mean.vec = matrix(0, ncol = dim), cov.mat = diag(dim), batch.size = 1000){
  # Input validation
  if (!is.numeric(dim) || length(dim) != 1 || dim < 1 || dim != floor(dim)) {
    stop(paste0("'dim' must be a positive integer, got ", deparse(dim)), call. = FALSE)
  }
  if (!is.numeric(mean.vec)) {
    stop(paste0("'mean.vec' must be numeric, got '", paste(class(mean.vec), collapse = "/"), "'"), call. = FALSE)
  }
  mean.len <- if (is.matrix(mean.vec)) ncol(mean.vec) else length(mean.vec)
  if (mean.len != dim) {
    stop(paste0("'mean.vec' length must equal 'dim' (expected ", dim, ", got ", mean.len, ")"), call. = FALSE)
  }
  if (!is.matrix(cov.mat) || !is.numeric(cov.mat)) {
    stop(paste0("'cov.mat' must be a numeric matrix, got '", paste(class(cov.mat), collapse = "/"), "'"), call. = FALSE)
  }
  if (nrow(cov.mat) != dim || ncol(cov.mat) != dim) {
    stop(paste0("'cov.mat' must be a ", dim, "x", dim, " matrix, got ", nrow(cov.mat), "x", ncol(cov.mat)), call. = FALSE)
  }

  mean.vec <- as.numeric(mean.vec)
  batch <- NULL
  batch.idx <- 0L
  batch.adjuster <- NULL

  addedHighDimGaussian <- function(point, covariance.adjuster = 1, dim = ""){
    if (is.null(batch) || batch.idx >= nrow(batch) || !identical(covariance.adjuster, batch.adjuster)) {
      batch <<- mvtnorm::rmvnorm(batch.size, mean = mean.vec, sigma = covariance.adjuster * cov.mat)
      batch.idx <<- 0L
      batch.adjuster <<- covariance.adjuster
    }
    batch.idx <<- batch.idx + 1L
    point + batch[batch.idx, ]
  }

}
