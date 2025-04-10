mclustDensityFunction <- function(env.model = NULL, presence.model = NULL, dim = "", threshold = 0.01){

  densityAtPointEstimator <- function(point){
    point <- as.matrix(sf::st_drop_geometry(point[dim]))
    density <- mclust::predict.densityMclust(env.model, point)
    if (density < threshold ) return(0)
    return(max(0, 1 - mclust::predict.densityMclust(presence.model, point)  * 5)) # this approach is dubeaous
    }
  return(densityAtPointEstimator)
}
