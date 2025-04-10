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
mcmcSampling <- function(dataset = NULL, dimensions= list(""), densityFunction = alwaysOne,
                         proposalFunction = addHighDimGaussian(dim = length(dimensions)), n.sample.points = 0){
  pb <- utils::txtProgressBar(min = 0, max = n.sample.points, style = 3)
  starting.index <- stats::runif(1,1,nrow(dataset))
  current.point <- dataset[starting.index,]
  dataset <- dataset[-starting.index,]
  sampled.points <- dataset[0, ]
  points.rejected <- 0
  while (nrow(sampled.points) < n.sample.points) {

    proposed.point <- proposalFunction(current.point, dim = dimensions)
    if (acceptNextPoint(current.point, proposed.point, densityFunction)){
      current.point <- proposed.point
      utils::setTxtProgressBar(pb, nrow(sampled.points))
      sampled.points <- rbind(sampled.points, proposed.point)
    }
    else points.rejected <- points.rejected + 1
    cat("\rPoints rejected:", points.rejected, "Points accepted:", nrow(sampled.points))
  }
  cat("Points rejected: ", points.rejected)
  #sampled.points <- apply(sampled.points, 1, function(point) mapBackOnRealPoints(dataset, point, dimensions))
  return(sampled.points)
}
