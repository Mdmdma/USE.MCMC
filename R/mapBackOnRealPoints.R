#' Map Back on Real Points
#' searches the closest point in the dataset regarding the given point and the dimensions given
#'
#' @param dataset dataframe of the target dataset
#' @param point dataframe containing the given point
#' @param dim string vector containing the names of the dimensions that should be respected
#' @param threshold threshold obove which we consider the distance too big and want to discard the point
#'
#' @returns closest point that we could find in the dataset with the distance added as a parameter
#' @export
#'
mapBackOnRealPoints <- function(dataset = NULL, point = NULL, dim ="", threshold = 0.5){
  distances <- apply(dataset, 1, function(row) euclidianMetric(row, point, dim))
  min.dist.index <- which.min(distances)
  min.dist <- distances[min.dist.index]
  best.point <- dataset[min.dist.index,]
  best.point$distanceFromSampledPoint <- min.dist
  if (min.dist < threshold) return(best.point)
  return(NA)
}
