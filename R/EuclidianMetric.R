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
euclidianMetric <- function(point.a=NULL, point.b =NULL, dim = ""){
  point.a <- sf::st_drop_geometry(point.a)
  point.b <- sf::st_drop_geometry(point.b)
  m <- Map(`-`, point.a[dim], point.b[dim]) %>% unlist() %>% abs() %>% sum()
}
