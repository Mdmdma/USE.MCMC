plotDensity2dpro <- function(dataset, xlim = c(0,1), ylim = c(0,1),
                          densityFunction = NULL, resolution = 10,
                          colorPalette = "viridis") {
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
  if (colorPalette == "viridis") {
    if (!requireNamespace("viridisLite", quietly = TRUE)) {
      colors <- heat.colors(100)
    } else {
      colors <- viridisLite::viridis(100)
    }
  } else if (colorPalette == "heat") {
    colors <- heat.colors(100)
  } else if (colorPalette == "topo") {
    colors <- topo.colors(100)
  } else {
    # Default custom gradient
    colors <- colorRampPalette(c("white", "yellow", "orange", "red"))(100)
  }

  # Plot with image
  image(x = x_seq, y = y_seq, z = density_matrix,
        col = colors,
        xlab = "PC1", ylab = "PC2",
        main = "Density Plot",
        xlim = xlim, ylim = ylim)

  # Add contour lines for better visualization
  contour(x = x_seq, y = y_seq, z = density_matrix,
          add = TRUE, col = "black", lwd = 0.5, alpha = 0.3)

  # Return the grid data invisibly (useful for further analysis)
  #points(dataset$PC1, dataset$PC2, pch = 20)
  invisible(list(grid = grid, matrix = density_matrix))
}

