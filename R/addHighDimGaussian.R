#' addHighDimGaussian
#'
#' Adds high dimensional Gaussian noise to specific parameters of a point given as a dataframe.
#'
#'
#' @param dim integer that specifies the number of dimensions that should be altered
#' @param mean_vec vector of the means of the Gaussian to be added
#' @param cov_mat covariance matrix of the Gaussian to be added
#'
#' @returns a function that takes a point given as a dataframe as input and returns it with gaussian noise added to the specified dimensions
#' @export
#' @importFrom dplyr %>%
#'
addHighDimGaussian <- function(dim = 0, mean_vec = matrix(0, ncol = dim), cov_mat = diag(dim)){
  addedHighDimGaussian <- function(currentPoint, dim = ""){
    randomVector <- mvtnorm::rmvnorm(1, mean = mean_vec, sigma = cov_mat) %>% c(0)
    currentPoint <- currentPoint %>% dplyr::mutate(dplyr::across(dplyr::all_of(dim), ~ . + randomVector[match(dplyr::cur_column(), dim)]))
  }

}
