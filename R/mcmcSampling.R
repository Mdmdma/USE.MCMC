mcmcSampling <- function(dataset = NULL, sampling.extend = list(c(0,1)), desityFunction = alwaysOne, proposalFunction = highDimGaussian(dim=length(sampling.extend)) , n.sample.points = 0, dimensions= list("") ){
  # Check parameters to ensure correct input
  if (!is.list(sampling.extend)) {
    stop("'sampling.extend' must be a list.")
  }
  if (!is.list(dimensions)) {
    stop("'dimensions' must be a list.")
  }
  if (length(sampling.extend) != length(dimensions)) {
    stop("'sampling.extend' and 'dimensions' must be lists of the same size.")
  }



}
