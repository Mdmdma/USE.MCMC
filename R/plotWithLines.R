plot_points_with_lines <- function(df, cols) {
  if (length(cols) < 2) {
    stop("At least two columns are needed for plotting.")
  }

  # Extract x and y coordinates
  x <- df[[cols[1]]]
  y <- df[[cols[2]]]

  # Plot points
  plot(x, y, type = "b", pch = 16, col = grDevices::hsv(seq_along(x)/length(x), 1,1, alpha = 1),
       xlab = cols[1], ylab = cols[2], main = "Connected Data Points")

  # Add labels
  graphics::text(x, y, labels = seq_along(x), pos = 3, cex = 0.8)
}
