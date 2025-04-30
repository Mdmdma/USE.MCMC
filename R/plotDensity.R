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
plotDensity <- function(dataset, species = NULL, xlim = c(0,1), ylim = c(0,1),
                        densityFunction = NULL, resolution = 10) {
  if(is.null(dataset)) stop("at least one line of data has to be supplied as a template for the density function")

  # Create a grid matrix
  x_seq <- seq(xlim[1], xlim[2], length.out = resolution)
  y_seq <- seq(ylim[1], ylim[2], length.out = resolution)
  grid <- expand.grid(PC1 = x_seq, PC2 = y_seq)

  # Create a template from sample.point with correct columns
  sample.point <- dataset[1,]
  sample.point$density <- 0
  template <- sample.point[1,]

  # Pre-allocate density values
  density_values <- numeric(nrow(grid))

  # Calculate densities
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

  # Create a matrix for contour plotting
  density_matrix <- matrix(grid$density, nrow = length(x_seq), ncol = length(y_seq))
  colnames(density_matrix) <- y_seq
  rownames(density_matrix) <- x_seq

  # We'll just use the grid directly, no need to reshape with reshape2

  # Create base plot using the grid data directly
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = PC1, y = PC2, z = density))

  # Add filled contours with viridis colors
  p <- p + ggplot2::geom_tile(ggplot2::aes(fill = density))
  p <- p + ggplot2::scale_fill_viridis_c()

  # Add contour lines
  p <- p + ggplot2::geom_contour(color = "black", alpha = 0.3, linewidth = 0.5)

  # Add points from dataset if available
  if (nrow(dataset) > 1) {
    p <- p + ggplot2::geom_point(data = dataset,
                                 ggplot2::aes(x = PC1, y = PC2, color = "Dataset"),
                                 inherit.aes = FALSE)
  }

  # Add species points in red if provided
  if (!is.null(species)) {
    p <- p + ggplot2::geom_point(data = species,
                                 ggplot2::aes(x = PC1, y = PC2, color = "Species"),
                                 inherit.aes = FALSE)
  }

  # Add color scale for points
  p <- p + ggplot2::scale_color_manual(name = "Points",
                                       values = c("Dataset" = "black", "Species" = "red"))

  # Set plot limits and labels
  p <- p + ggplot2::xlim(xlim) +
    ggplot2::ylim(ylim) +
    ggplot2::labs(x = "PC1", y = "PC2", title = "Density Plot",
                  fill = "Density")  # Label for the fill legend

  # Add theme elements
  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # Print the plot
  print(p)

  # Return the grid data invisibly (useful for further analysis)
  invisible(list(grid = grid, matrix = density_matrix, plot = p))
}
