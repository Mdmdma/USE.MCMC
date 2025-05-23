#' MCMC sampling from a given dataset
#'
#' @param dataset sf dataframe from which the points are sampled
#' @param dimensions string vector containing the dimensions that should be included in the random walk
#' @param densityFunction Function that can take a point given as a sf dataframe as a input and returns the target density at that location.
#' @param proposalFunction Function that can take a point given as a sf dataframe and a vector of strings specifying the row names that should be changed as a input and returns a new proposed point
#' @param n.sample.points Number of points to be sampled
#' @param burnIn Boolean sets if a burn in shoud be done before the sampling
#' @param verbose Boolean to toggle progress updates
#' @returns A sf dataframe containing the sampled points
#' @export
#'
mcmcSampling <- function(dataset = NULL, dimensions= list(""), densityFunction = alwaysOne,
                         proposalFunction = addHighDimGaussian(dim = length(dimensions)), n.sample.points = 0, burnIn = TRUE, verbose = TRUE){

  starting.index <- stats::runif(1, 1, nrow(dataset))
  current.point <- dataset[starting.index, ]
  current.point <- sf::st_drop_geometry(current.point)
  current.point$density <- densityFunction(current.point)
  # burn in
  covariance.correction <- 1
  if(burnIn) {
    cat("Burn in\n")
    points.rejected <- 0
    num.burnin.samples <- 1000
    if (verbose){
      pb.burnin <- utils::txtProgressBar(min = 0, max = num.burnin.samples, style = 3)
    }
   while (points.rejected / num.burnin.samples < 0.21| points.rejected / num.burnin.samples > 0.25) {
      # the numbers of the condition depend on the exact threshold, Gelman, Roberts, and Gilks (1996) proposes 0.23 was optimal
      points.rejected <- 0
      for (i in 1:num.burnin.samples) {
        proposed.point <- proposalFunction(current.point, covariance.adjuster = covariance.correction, dim = dimensions)
        proposed.point$density <- densityFunction(proposed.point)
        if (acceptNextPoint(current.point, proposed.point)){
          current.point <- proposed.point
        }
        else points.rejected <- points.rejected + 1

        if (verbose){
          cat("\rCurrent rejection ratio: ", points.rejected / i)
          utils::setTxtProgressBar(pb.burnin, i)
        }
      }
      if (points.rejected / num.burnin.samples < 0.21) covariance.correction <- covariance.correction * stats::rnorm(1, mean = 1.3, sd = 0.1)
      if (points.rejected /num.burnin.samples > 0.25) covariance.correction <- covariance.correction * stats::rnorm(1, mean = 0.7, sd = 0.1)
      if (verbose){
        cat("\rThe current rejection rate is", points.rejected /num.burnin.samples, ", the currend covariance adjustment factor is ", covariance.correction, "\n")
      }

    }
  }
  else cat("Burn in skipped")

  cat("\nThe final covariance correction is ", covariance.correction, "\n")
  sampled.points <- data.frame(matrix(NA, nrow = n.sample.points, ncol = ncol(current.point)))
  colnames(sampled.points) <- colnames(current.point)
  points.rejected <- 0
  points.accepted <- 0

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n.sample.points, style = 3)
  }

  while (points.accepted < n.sample.points) {
    proposed.point <- proposalFunction(current.point, covariance.adjuster = covariance.correction, dim = dimensions)
    proposed.point$density <- densityFunction(proposed.point)
    points.accepted <- points.accepted + 1
    if (acceptNextPoint(current.point, proposed.point)){
      current.point <- proposed.point
      sampled.points[points.accepted, ] <- current.point
    }
    else {
      sampled.points[points.accepted, ] <- current.point
      points.rejected <- points.rejected + 1
    }
    if (verbose){
      cat("\rPoints rejected:", points.rejected, "Points accepted:", points.accepted)
      utils::setTxtProgressBar(pb, points.accepted)
    }
  }
  cat("\nPoints rejected: ", points.rejected)

  return(sampled.points)
}

