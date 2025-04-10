#' Function to plot the density of a density
#'
#'
#' @param dataset dataframe containing points that can be added to the desity field. The parameter is needed as the code
#' uses it as a template for points. If its size is greater than one the points are added to the plot
#'
#' @param species dataframe containing the species points to be added to the plot
#' @param xlim Limits of the plot
#' @param ylim Limits of the plot
#' @param densityFunction Density to be plottet. The function should take a point as its argument and give back a float
#' @param resolution number of gridcells in each direction
#'
#' @returns list containing the sampling grid as well as the density matrix
#' @export
plotDensity2dpro <- function(dataset, species = NULL, xlim = c(0,1), ylim = c(0,1),
                          densityFunction = NULL, resolution = 10) {
  # Create a grid matrix instead of nested loops
  x_seq <- seq(xlim[1], xlim[2], length.out = resolution)
  y_seq <- seq(ylim[1], ylim[2], length.out = resolution)
  grid <- expand.grid(PC1 = x_seq, PC2 = y_seq)

  # Create a template from sample.point with correct columns
  sample.point <- dataset[1,]
  sample.point$density <- 0
  template <- sample.point[1,]

  # Pre-allocate density values
  density_values <- numeric(nrow(grid))

  # Calculate densities (vectorized if possible)
  if (is.function(densityFunction)) {
    density_values <- mapply(function(pc1, pc2) {
      template$PC1 <- pc1
      template$PC2 <- pc2
      densityFunction(template)
    }, grid$PC1, grid$PC2)
  } else {
    stop("A density function must be provided")
  }

  # Add density to grid
  grid$density <- density_values

  # Create a matrix for image plotting (much faster than points)
  density_matrix <- matrix(grid$density, nrow = length(x_seq), ncol = length(y_seq))

  # Plot using image() which is much faster and smoother than points
  # Choose color palette

  colors <- grDevices::gray.colors(100)

  # Plot with image
  graphics::image(x = x_seq, y = y_seq, z = density_matrix,
        col = colors,
        xlab = "PC1", ylab = "PC2",
        main = "Density Plot",
        xlim = xlim, ylim = ylim)

  # Add contour lines for better visualization
  graphics::contour(x = x_seq, y = y_seq, z = density_matrix,
          add = TRUE, col = "black", lwd = 0.5, alpha = 0.3)

  # Return the grid data invisibly (useful for further analysis)
  if (nrow(dataset)>1) graphics::points(dataset$PC1, dataset$PC2, pch = 20)
  if (!is.null(species)) graphics::points(species$PC1, species$PC2, pch = 20, col = "red")
  invisible(list(grid = grid, matrix = density_matrix))
}

