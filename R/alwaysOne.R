#' AlwaysOne
#'
#' As the name implies, this function returns always 1, no matter the input. It can be used as a dummy desity function to perform uniform sampling of higher dimensional space. As it works with any Input, it can also be used as a placeholder to evade errors.
#'
#' @param ... Any input is possible, as it does not have any effect
#'
#' @returns 1, always
#' @keywords internal
#' alwaysOne(stop("even this input is fine"))
alwaysOne <- function(...){
  1
}
