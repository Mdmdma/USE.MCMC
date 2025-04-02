#' Plott a dequence of points in 2d
#'
#' @param df dataframe containing the observation
#' @param cols vector containing the names of the columns to be plotted
#' @param limits sets the limits of the plot. This is usefull to ensure comparable plot sizes
#' @param title title of the plot
#'
#' @returns the desired plot
#' @export
#'
plot_points_with_lines <- function(df, cols, limits = NULL, title = "Connected Data Points") {
  if (length(cols) < 2) {
    stop("At least two columns are needed for plotting.")
  }

  # Extract x and y coordinates
  x <- df[[cols[1]]]
  y <- df[[cols[2]]]

  # Generate colors for each point
  point_colors <- grDevices::hsv(seq_along(x) / length(x), 1, 1, alpha = 0.5)

  # Determine axis limits
  if (!is.null(limits) && length(limits) == 2) {
    xlim <- limits[[1]]
    ylim <- limits[[2]]
  } else {
    xlim <- range(x, na.rm = TRUE)
    ylim <- range(y, na.rm = TRUE)
  }

  # Plot empty canvas with specified limits
  plot(x, y, type = "n", xlab = cols[1], ylab = cols[2],
       main = title, xlim = xlim, ylim = ylim)

  # Draw the lines with the same colors as the points
  for (i in seq_along(x)[-1]) {
    segments(x[i-1], y[i-1], x[i], y[i], col = point_colors[i], lwd = 2)
  }

  # Plot the points
  points(x, y, pch = 16, col = point_colors)

  # Add labels
  graphics::text(x, y, labels = seq_along(x), pos = 3, cex = 0.8)
}

