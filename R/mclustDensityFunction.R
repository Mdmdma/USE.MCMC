mclustDensityFunction <- function(model, dim = "", threshold = 0.01){

  densityAtPointEstimator <- function(point){
    point <- as.matrix(sf::st_drop_geometry(point[dim]))
    density <- mclust::predict.densityMclust(model, point)
    if (density < threshold ) return(0)
    return(1)
  }
  return(densityAtPointEstimator)
}
