#' Helper to create a Density function that uses mclust Gaussian mixtures
#'
#' @param env.model mclust gaussian mixture that uses points
#' @param presence.model mclust gaussian mixture that uses points
#' @param dim string vector specifing the names of the dimensions
#' @param threshold sets the curoff density from the environment
#'
#' @returns Function that can calculates the density at a point
#' @export
#'
mclustDensityFunction <- function(env.model = NULL, species.model = NULL, dim = "", threshold = 0.01, species.cutoff.threshold = 0.1){

  densityAtPointEstimator <- function(point){
    point <- as.matrix(sf::st_drop_geometry(point[dim]))
    density <- mclust::predict.densityMclust(env.model, point)
    if (density < threshold ) return(threshold / 1000)
    return(max(threshold / 1000, 1 - mclust::predict.densityMclust(species.model, point)  / species.cutoff.threshold)) # this approach is dubeaous
    }
  return(densityAtPointEstimator)
}
