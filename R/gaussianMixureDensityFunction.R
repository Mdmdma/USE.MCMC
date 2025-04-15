gaussianMixtureDensityFunction <- function(env.data, dim = "", threshold = 0.01){

  env.data <- sf::st_drop_geometry(env.data[dim])
  print("Fitting environmental conditions")
  environmental.data.model <- mclust::densityMclust(env.data, plot = FALSE)
  mclust::plot.densityMclust(environmental.data.model, what = "density")
  summary(environmental.data.model)

  densityAtPointEstimator <- function(point){
    point <- as.matrix(sf::st_drop_geometry(point[dim]))
    density <- mclust::predict.densityMclust(environmental.data.model, point)
    if (density < threshold ) return(threshold/ 1000)
    return(1)
  }
  return(densityAtPointEstimator)
}
