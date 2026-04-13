#' Accept next point ?
#'
#' \code{acceptNextPoint} checks whether the proposed point should be accepted
#'
#' @param current.density Numeric scalar, density at the current point.
#' @param proposed.density Numeric scalar, density at the proposed point.
#'
#' @returns TRUE if the proposed point should be accepted
#' @keywords internal
#'
acceptNextPoint <- function(current.density, proposed.density){
  acceptance.ratio <- proposed.density / current.density
  if (is.nan(acceptance.ratio) || is.na(acceptance.ratio)) return(FALSE)
  if (acceptance.ratio > 1) return(TRUE)
  if (stats::runif(1) < acceptance.ratio) return(TRUE)
  return(FALSE)
}
