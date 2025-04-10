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
plotPointsWithLines <- function(df, cols, limits = NULL, title = "Connected Data Points") {
  if (length(cols) < 2) {
    stop("At least two columns are needed for plotting.")
  }

  # Extract x and y coordinates
  x <- df[[cols[1]]]
  y <- df[[cols[2]]]
  v = 1
  if (is.null(df[[cols[4]]])){
    if (min(df[[cols[3]]]) != max(df[[cols[3]]])){
      v <- 1-(df[[cols[3]]]-min(df[[cols[3]]]))/(max(df[[cols[3]]])-min(df[[cols[3]]]))
    }
  }

  # Generate colors for each point
  point.colors <- grDevices::hsv(seq_along(x) / length(x), 1, v, alpha = 0.2)

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
    graphics::segments(x[i-1], y[i-1], x[i], y[i], col = point.colors[i], lwd = 2)
  }

  # Plot the points
  graphics::points(x, y, pch = 16, col = point.colors)

  # Add labels
  #graphics::text(x, y, labels = seq_along(x), pos = 3, cex = 0.8)
}
