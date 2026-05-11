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
#' @returns maximal distance a point should be remapped from to have originated from a region inside of the environmental space.
#' @export
optimalDistanceThresholdNn <- function(env.data = NULL,
                                      index.for.cutof = 5,
                                      dimensions = c("PC1", "PC2"),
                                      num.neighbors = 3){
  # Input validation
  if (is.null(env.data)) {
    stop("'env.data' must be provided (got NULL)", call. = FALSE)
  }
  if (!is.data.frame(env.data) && !inherits(env.data, "sf")) {
    stop(paste0("'env.data' must be a data.frame or sf object, got '",
                paste(class(env.data), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.character(dimensions) || length(dimensions) < 1) {
    stop("'dimensions' must be a character vector of column names", call. = FALSE)
  }
  env_cols <- if (inherits(env.data, "sf")) setdiff(names(env.data), attr(env.data, "sf_column")) else names(env.data)
  missing_dims <- setdiff(dimensions, env_cols)
  if (length(missing_dims) > 0) {
    stop(paste0("'dimensions' contains columns not found in 'env.data': ",
                paste0("'", missing_dims, "'", collapse = ", "),
                ". Available columns: ", paste(env_cols, collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(index.for.cutof) || length(index.for.cutof) != 1 || index.for.cutof < 1 || index.for.cutof != floor(index.for.cutof)) {
    stop(paste0("'index.for.cutof' must be a positive integer, got ", deparse(index.for.cutof)), call. = FALSE)
  }
  if (!is.numeric(num.neighbors) || length(num.neighbors) != 1 || num.neighbors < 1 || num.neighbors != floor(num.neighbors)) {
    stop(paste0("'num.neighbors' must be a positive integer, got ", deparse(num.neighbors)), call. = FALSE)
  }

  env.data.cleaned <- sf::st_drop_geometry(env.data)
  nearest.neighbors.distance <- FNN::knn.dist(env.data.cleaned[dimensions],
                                              k = num.neighbors) %>%
    as.vector()
  # Partial sort: we only need the `index.for.cutof`-th largest distance, so avoid
  # a full O(n log n) sort over what may be a length-(n * num.neighbors) vector.
  top.cutoff <- -sort(-nearest.neighbors.distance, partial = index.for.cutof)[index.for.cutof]
  distance.threshold <- top.cutoff / 2
  return(distance.threshold)
}

