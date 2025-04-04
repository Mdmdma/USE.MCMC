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
mcmcSampling <- function(dataset = NULL, dimensions= list(""), densityFunction = alwaysOne, proposalFunction = addHighDimGaussian(dim = lengt(dimensions)), n.sample.points = 0){
  pb <- txtProgressBar(min = 0, max = n.sample.points, style = 3)
  starting.index <- stats::runif(1,1,nrow(dataset))
  current.point <- dataset["25918",]
  dataset <- dataset[-starting.index,]
  sampled.points <- dataset[0, ]
  points.rejected <- 0
  while (nrow(sampled.points) < n.sample.points) {

    proposed.point <- proposalFunction(current.point, dim = dimensions)
    if (acceptNextPoint(current.point, proposed.point, densityFunction)){
      #distances <- apply(dataset, 1, function(row) euclidianMetric(row, proposed.point, dimensions))
      #min.dist.index <- which.min(distances)
      #current.point <- dataset[min.dist.index,]
      #sampled.points <- rbind(sampled.points, current.point)
      current.point <- proposed.point
      setTxtProgressBar(pb, nrow(sampled.points))
      sampled.points <- rbind(sampled.points, proposed.point)
    }
    else points.rejected <- points.rejected + 1
    cat("\rPoints rejected:", points.rejected, "Points accepted:", nrow(sampled.points))
  }
  cat("Points rejected: ", points.rejected)
  return(sampled.points)
}
