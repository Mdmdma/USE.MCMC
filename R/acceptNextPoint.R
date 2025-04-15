#' Accept next point ?
#'
#' \code{acceptNextPoint} checks whether the proposed point should be accepted
#'
#' @param current.point Current state of the markov chain, given as a sf dataframe with a column called $density containing the value of the target density.
#' @param proposed.point Proposed next state of the markov chain, given as a sf dataframe with a column called $density containing the value of the target density.
#'
#' @returns TRUE if the proposed point should be accepted
#' @keywords internal
#'
acceptNextPoint <- function(current.point, proposed.point){
  acceptance.ratio <- proposed.point$density / current.point$density
  if (acceptance.ratio > 1) return(TRUE)
  if (stats::runif(1) < acceptance.ratio) return(TRUE)
  return(FALSE)
}
