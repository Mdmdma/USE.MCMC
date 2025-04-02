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
#'
euclidianMetric <- function(pointA=NULL, pointB =NULL, dim = ""){
  pointA <- sf::st_drop_geometry(pointA)
  pointB <- sf::st_drop_geometry(pointB)
  m <- Map(`-`, pointA[dim], pointB[dim]) %>% unlist() %>% sum() %>% abs()
}
