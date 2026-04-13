gaussianMixtureDensityFunction <- function(env.data, dim = "", threshold = 0.01){

  env.data <- sf::st_drop_geometry(env.data[dim])
  print("Fitting environmental conditions")
  environmental.data.model <- mclust::densityMclust(env.data, plot = FALSE)
  mclust::plot.densityMclust(environmental.data.model, what = "density")
  summary(environmental.data.model)

  # Precompute GMM parameters for fast evaluation
  env.pre <- precompute_gmm_params(environmental.data.model)
  threshold.floor <- threshold / 1000

  densityAtPointEstimator <- function(point){
    density <- fast_gmm_density(point, env.pre)
    if (density < threshold) return(threshold.floor)
    return(1)
  }
  return(densityAtPointEstimator)
}
