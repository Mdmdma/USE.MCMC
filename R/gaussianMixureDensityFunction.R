gaussianMixtureDensityFunction <- function(environmentalData, dim = "", threshold = 0.01){
  environmentalData <- sf::st_drop_geometry(environmentalData[dim])
  print("Fitting environmental conditions")
  environmental.data.model <- mclust::densityMclust(environmentalData, plot = FALSE)
  mclust::plot.densityMclust(environmental.data.model, what = "density")
  summary(environmental.data.model)

  densityAtPointEstimator <- function(point){
    point <- as.matrix(sf::st_drop_geometry(point[dim]))
    density <- mclust::predict.densityMclust(environmental.data.model, point)
    if (density < threshold ) return(0)
    return(1)
  }
  return(densityAtPointEstimator)
}
