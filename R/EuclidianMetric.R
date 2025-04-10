#' Euclidean metric on specified columns
#'
#' \code{euclidianMetric} calculates the euclidean metric of specified columns of two points given as dataframes
#'
#' @param pointA First point
#' @param pointB Second point
#' @param dim vector containing the names of the columns that should be included in the computation
#'
#' @returns euclidean metric of the two points
#' @keywords internal
#' @importFrom dplyr %>%
euclidianMetric <- function(pointA=NULL, pointB =NULL, dim = ""){
  point.A <- sf::st_drop_geometry(pointA)
  point.B <- sf::st_drop_geometry(pointB)
  m <- Map(`-`, point.A[dim], point.B[dim]) %>% unlist() %>% abs() %>% sum()
}
