#' Accept next point ?
#'
#' \code{acceptNextPoint} checks whether the proposed point should be accepted
#'
#' @param current.point Current state of the markov chain, given as a sf dataframe.
#' @param proposed.point Proposed next state of the markov chain, given as a sf dataframe.
#' @param densityFunction Function that can take a point given as a sf dataframe as a input and returns the target density at that location.
#'
#' @returns TRUE if the proposed point should be accepted
#' @keywords internal
#'
acceptNextPoint <- function(current.point, proposed.point, densityFunction){
  acceptance.ratio <- densityFunction(proposed.point) / densityFunction(current.point)
  if (stats::runif(1) < acceptance.ratio) return(TRUE)
  return(FALSE)
}
