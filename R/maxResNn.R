#' maxResNN
#' is a function that can be used to compute a reasonable grid resolution for nearest neighbor based uniform sampling.
#' The core idea behind its working principle is that we want to expect a grid cell to contain points if it overlaps with the environment.
#' This implementation looks at the low density regions using distance to n neighbors as a proxy.
#' The approach assumes that both coordinates have similar range, as the axes are not weighted when computing Neighbors and converting from distances to number of grid cells
#' If multiple points from the same grid cell should be sampled, the number of neighbors included in the computation should be set accordingly
#'
#' @param env.data.raster Raster containing the environmental parameters
#' @param dimensions vector containing the dimensions that should be used for the grid computation
#' @param low.end.of.inclueded.points low cutoff
#' @param high.end.of.included.points high cutoff
#' @param n.neighbors number of neighbors used for the computation
#' @param PCA can be set to true if rastPCA was already used to perform a pca to save time recomputing
#'
#' @returns maximal number of grid cells useful for the nearest neighbor based approach
#' @export
#'
#'
maxResNn<- function(env.data.raster, dimensions = c("PC1", "PC2") , low.end.of.inclueded.points = 20, high.end.of.included.points = 4 , n.neighbors = 2, PCA = FALSE) {
  # Input validation
  if (!is.logical(PCA) || length(PCA) != 1) {
    stop(paste0("'PCA' must be a single logical value, got '", paste(class(PCA), collapse = "/"), "'"), call. = FALSE)
  }
  if (!PCA) {
    check_raster_input(env.data.raster, "env.data.raster")
  } else {
    if (!is.list(env.data.raster) || is.null(env.data.raster$PCs)) {
      stop("When PCA=TRUE, 'env.data.raster' must be a list with a '$PCs' element (result of rastPCA)", call. = FALSE)
    }
  }
  if (!is.character(dimensions) || length(dimensions) < 2) {
    stop("'dimensions' must be a character vector with at least 2 elements", call. = FALSE)
  }
  if (!is.numeric(low.end.of.inclueded.points) || length(low.end.of.inclueded.points) != 1 || low.end.of.inclueded.points < 1) {
    stop(paste0("'low.end.of.inclueded.points' must be a positive number, got ", deparse(low.end.of.inclueded.points)), call. = FALSE)
  }
  if (!is.numeric(high.end.of.included.points) || length(high.end.of.included.points) != 1 || high.end.of.included.points < 1) {
    stop(paste0("'high.end.of.included.points' must be a positive number, got ", deparse(high.end.of.included.points)), call. = FALSE)
  }
  if (!is.numeric(n.neighbors) || length(n.neighbors) != 1 || n.neighbors < 1 || n.neighbors != floor(n.neighbors)) {
    stop(paste0("'n.neighbors' must be a positive integer, got ", deparse(n.neighbors)), call. = FALSE)
  }

  # PCA if a pca has been performed upstream
  if(PCA) {
    rpc <- env.data.raster
  } else {
    rpc <- rastPCA(env.data.raster, stand = TRUE)
  }


  # convert to dataframe
  pc.df <- rpc$PCs %>%
    as.data.frame(xy = TRUE) %>%
    na.omit() %>%
    dplyr::select(dimensions)

  neighbor.data <- FNN::get.knnx(pc.df[dimensions], pc.df[dimensions], k = 1 + n.neighbors) # plus one as each point is its own closest <<neighbor>>
  neighbor.index <- neighbor.data$nn.index
  # Columns 2:(n.neighbors+1) are the distances to the n.neighbors nearest non-self
  # points (column 1 is the self-distance). The previous `2:n.neighbors + 1` was
  # parsed as `(2:n.neighbors) + 1` and silently dropped the nearest non-self.
  neighbor.distance <- neighbor.data$nn.dist[, 2:(n.neighbors + 1), drop = FALSE]
  neighbor.distance.flat <- as.vector(neighbor.distance)
  sorted.distances <- sort(neighbor.distance.flat)
  top.end.points.without.outlayers.distance <- sorted.distances[(length(sorted.distances) - low.end.of.inclueded.points):(length(sorted.distances) - high.end.of.included.points)]
  data.ranges <- sapply(pc.df, function(col) range(col))
  max.dim.range <- max(data.ranges[2, 1] - data.ranges[1, 1],
                       data.ranges[2, 2] - data.ranges[1, 2])
  max.num.of.cells <- ceiling(max.dim.range / mean(top.end.points.without.outlayers.distance))

  return(max.num.of.cells)
}
