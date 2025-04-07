mclustDensityFunction <- function(environmentalModel = NULL, presenceModel = NULL, dim = "", threshold = 0.01){

  densityAtPointEstimator <- function(point){
    point <- as.matrix(sf::st_drop_geometry(point[dim]))
    density <- mclust::predict.densityMclust(environmentalModel, point)
    if (density < threshold ) return(0)
    return(1 - mclust::predict.densityMclust(presenceModel, point)  * 5)
    }

  return(densityAtPointEstimator)
}
