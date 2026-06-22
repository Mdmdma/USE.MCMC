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
#'  \strong{Dimension dependence.} The bare \code{D_low / 2} rule is only correct in
#'  two dimensions. The deepest \emph{legitimate} interior gap of a point cloud (the
#'  Voronoi-cell circumradius / largest empty ball that still sits inside the
#'  support) grows with dimension relative to the local point spacing: for a cloud
#'  of spacing \eqn{a} it scales like \eqn{(\sqrt{d}/2)\,a}. Keeping the fixed
#'  \code{/2} factor as \code{d} grows therefore \emph{over-rejects} genuine
#'  low-density interior points and punches spurious holes in exactly the sparse
#'  region uniform sampling is meant to fill. \code{dim.correction} rescales the
#'  threshold by a factor that is \eqn{1} at \code{d = 2} (so two-dimensional
#'  behaviour is unchanged) and tracks this growth for \code{d > 2}. The NN
#'  \emph{remapping} itself needs no such correction — it is uniform over the
#'  support in every dimension.
#'
#' @param env.data Dataframe containing the environmental observation
#' @param index.for.cutof Index that is supposed to describe a low density point but not an outlayer
#' @param dimensions Vector containing the dimensions included in the analysis
#' @param num.neighbors Number of Neighbors
#' @param dim.correction Dimension correction for the rejection radius. Either a
#'   single positive numeric used as a literal multiplier, or one of
#'   \code{"voronoi"} (default; multiply by \eqn{\sqrt{d/2}}, the cubic deep-hole
#'   growth), \code{"simplex"} (multiply by \eqn{\sqrt{3d/(2(d+1))}}, a milder
#'   variant that saturates near \eqn{1.22} and is preferable if the default leaks
#'   pseudo-absences outside the support), or \code{"none"} (no correction, the
#'   legacy behaviour). All options equal \eqn{1} at \code{d = 2}.
#'
#' @returns maximal distance a point should be remapped from to have originated from a region inside of the environmental space.
#' @export
optimalDistanceThresholdNn <- function(env.data = NULL,
                                      index.for.cutof = 5,
                                      dimensions = c("PC1", "PC2"),
                                      num.neighbors = 3,
                                      dim.correction = "voronoi"){
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
  # Half the low-density spacing (2D-calibrated), rescaled by the dimension-aware
  # deep-hole factor so the threshold tracks the deepest legitimate interior gap as
  # the number of dimensions grows. The factor is 1 at d = 2.
  distance.threshold <- (top.cutoff / 2) * .nnThresholdDimFactor(length(dimensions), dim.correction)
  return(distance.threshold)
}

# Dimension correction factor for the NN distance threshold (internal). The deepest
# legitimate interior hole of a cloud with spacing a is the Voronoi-cell
# circumradius, ~ (sqrt(d)/2) * a for a cubic-lattice model; the legacy 2D rule used
# a/2, so to preserve the calibrated 2D behaviour we scale by sqrt(d/2) (== 1 at
# d = 2). NN-remapping is uniform-over-support in all d and needs no correction;
# only this support-membership threshold does. See optimalDistanceThresholdNn.
.nnThresholdDimFactor <- function(d, dim.correction = "voronoi") {
  if (is.numeric(dim.correction)) {
    if (length(dim.correction) != 1 || !is.finite(dim.correction) || dim.correction <= 0) {
      stop("numeric 'dim.correction' must be a single positive multiplier", call. = FALSE)
    }
    return(dim.correction)
  }
  dim.correction <- match.arg(dim.correction, c("voronoi", "simplex", "none"))
  switch(dim.correction,
         voronoi = sqrt(d / 2),                       # cubic deep-hole, ~sqrt(d)
         simplex = sqrt(3 * d / (2 * (d + 1))),       # regular-simplex, saturates ~1.22
         none    = 1)                                 # legacy: no correction
}

