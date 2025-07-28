#' Function to calculate the optimal distance threshold for nearest neighborhood sampling.
#'
#' The method is based on the assumption that the underlying dataset has a connected structure.
#'  By looking at the distance between points in areas with a low density, ie large distances between points, we can determine the distance that
#'  a realistic points could have if it originated from the observed environment to the realized data point. If we set the maximal distance to half the
#'  distance between points in the low density region, excluding outlayers, we exclude points that originate from "empty cells" in the classical case.
#'  We still over-sample the border regions, as the possible directions in which a possible real point can lay are limited compared to points on the inside.
#'
#'  A possible way to eliminate this issue could be to find n neighbors in the real dataset, look at the direction in which we can find them and reject
#'  a point with a probability of one minus the ration of the covered directions to the full hyper-sphere. In the high dimensional case this seems
#'  difficult, therefore it was not implemented.
#'
#'  The issue is analogue to  cells that only partially overlap with the environmental space. For the points in these cells the probability to be sampled
#'  is higher than it should be.
#'
#' @param env.data Dataframe containing the environmental observation
#' @param index.for.cutof Index that is supposed to describe a low density point but not an outlayer
#' @param dimensions Vector containing the dimensions included in the analysis
#' @param num.neighbors Number of Neighbors
#'
#' @returns
#' @export
#'
#' @examples
optimalDistanceThresholdNn <- function(env.data = NULL,
                                      index.for.cutof = 5,
                                      dimensions = c("PC1", "PC2"),
                                      num.neighbors = 3){

  env.data.cleaned <- sf::st_drop_geometry(env.data)
  nearest.neighbors.distance <- FNN::knn.dist(env.data.cleaned[dimensions],
                                              k = num.neighbors) %>%
    as.vector()
  sorted.nearest.neighbor.distances <- sort(nearest.neighbors.distance,
                                            decreasing=TRUE)
  distance.threshold <- sorted.nearest.neighbor.distances[2] / 2
  return(distance.threshold)
}

