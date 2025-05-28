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

  neighbor.data <- FNN::get.knnx(pc.df[dimensions], pc.df[dimensions],k = 1 + n.neighbors) # plus one as each points is its own closest <<neighbor>>
  neighbor.index <- neighbor.data$nn.index
  neighbor.distance <- neighbor.data$nn.dist[, 2:n.neighbors + 1] # remove the column with the distance to itself
  neighbor.distance.flat <- as.vector(neighbor.distance)
  sorted.distances <- sort(neighbor.distance.flat)
  top.end.points.without.outlayers.distance <- sorted.distances[(length(sorted.distances)- low.end.of.inclueded.points):(length(sorted.distances) - high.end.of.included.points)]
  data.ranges <- sapply(pc.df, function(col) range(col))
  data.ranges <- max(data.ranges[2,1] - data.ranges[1,1], data.ranges[2,2], data.ranges[2,1])
  max.num.of.cells <- ceiling(data.ranges / mean(top.end.points.without.outlayers.distance))

  return(max.num.of.cells)

  # interesting insights into distance related metrics
  graphics::par(mfrow = c(1,2))
  graphics::plot(sorted.distances)
  graphics::hist(sorted.distances)
}
