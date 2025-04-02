#' MCMC sampling from a given dataset
#'
#' @param dataset sf dataframe from which the points are sampled
#' @param dimensions string vector containing the dimensions that should be included in the random walk
#' @param densityFunction Function that can take a point given as a sf dataframe as a input and returns the target density at that location.
#' @param proposalFunction Function that can take a point given as a sf dataframe and a vector of strings specifying the row names that should be changed as a input and returns a new proposed point
#' @param n.sample.points Number of points to be sampled
#'
#' @returns A sf dataframe containing the sampled points
#' @export
#'
mcmcSampling <- function(dataset = NULL, dimensions= list(""), densityFunction = alwaysOne, proposalFunction = addHighDimGaussian , n.sample.points = 0){
  starting.index <- stats::runif(1,1,nrow(dataset))
  current.point <- dataset[starting.index,]
  dataset <- dataset[-starting.index,]
  sampled.points <- dataset[0, ]
  while (nrow(sampled.points) < n.sample.points) {
    print(paste("sampling point number", nrow(sampled.points)))
    proposed.point <- addHighDimGaussian(currentPoint = current.point, dim = dimensions)
    if (acceptNextPoint(current.point, proposed.point, densityFunction)){
      distances <- apply(dataset, 1, function(row) euclidianMetric(row, current.point, dimensions))
      min.dist.index <- which.min(distances)
      current.point <- dataset[min.dist.index,]
      dataset <- dataset[-min.dist.index,]
      sampled.points <- rbind(sampled.points, current.point)
      print(nrow(dataset))
    }
  }
  return(sampled.points)
}
