#' multivariate gaussian random number generator
#'
#' @param dim number of dimensions
#' @param mean_vec first moment
#' @param cov_mat second moment
#'
#' @returns a function that generates samples from a higher dimensional gaussian
#' @keywords internal

highDimGaussian <- function(dim = 1, mean_vec = matrix(1, ncol = dim), cov_mat = diag(dim)){
  rng <- function(){
    mvtnorm::rmvnorm(1, mean = mean_vec, sigma = cov_mat)
  }
}
