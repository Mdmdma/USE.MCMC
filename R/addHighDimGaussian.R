#' addHighDimGaussian
#'
#' Adds high dimensional Gaussian noise to specific parameters of a point given as a dataframe.
#'
#' @param currentPoint Current point as a dataframe from which a new point should be generated
#' @param dim vector with the names of the columns on which we want to add
#' @param mean_vec vector of the means of the Gaussian to be added
#' @param cov_mat covariance matrix of the Gaussian to be added
#'
#' @returns dataframe with Gaussian noise added to the specified columns
#' @keywords internal
#' @importFrom dplyr %>%
#'
addHighDimGaussian <- function(currentPoint, dim = "", mean_vec = matrix(0, ncol = length(dim)), cov_mat = 2.3 * diag(length(dim))){
  randomVector <- mvtnorm::rmvnorm(1, mean = mean_vec, sigma = cov_mat) %>% c(0)
  currentPoint <- currentPoint %>% dplyr::mutate(dplyr::across(dplyr::all_of(dim), ~ . + randomVector[match(cur_column(), dim)]))
}
