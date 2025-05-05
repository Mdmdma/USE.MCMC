#' Plott a dequence of points in 2d
#'
#' @param dataset dataframe containing the observation
#' @param cols vector containing the names of the columns to be plotted
#' @param limits sets the limits of the plot. This is usefull to ensure comparable plot sizes
#' @param title title of the plot
#'
#' @returns the desired plot
#' @export
#'
plotPointsWithLines <- function(dataset,
                                cols,
                                xlim = c(0,1), ylim = c(0,1),
                                title = "Connected Data Points") {
  p <- ggplot2::ggplot()
  if (length(cols) < 2) {
    stop("At least two columns are needed for plotting.")
  }

  # Create a copy of dataframe with row numbers for color gradient
  plot_df <- dataset
  plot_df$point_order <- 1:nrow(dataset)
  plot_df$point_order_normalized <- plot_df$point_order / nrow(dataset)

  # Initialize transparency value
  v <- 1
  # Check if we have a third column for transparency adjustment
  if (length(cols) >= 3 && !is.null(dataset[[cols[3]]]) && sum(is.na(dataset[[cols[3]]])) == 0) {
    if (min(dataset[[cols[3]]]) != max(dataset[[cols[3]]])) {
      # Same transparency calculation as original function
      plot_df$transparency <- 1 - (dataset[[cols[3]]] - min(dataset[[cols[3]]])) /
        (max(dataset[[cols[3]]]) - min(dataset[[cols[3]]])) / 2
      v <- plot_df$transparency
      # print(v)  # Keep the original debugging print
    }
  }

  # Create HSV colors for points and lines - with alpha
  hsv_colors_alpha <- character(nrow(plot_df))
  for(i in 1:nrow(plot_df)) {
    if(is.numeric(v) && length(v) > 1) {
      hsv_colors_alpha[i] <- grDevices::hsv(plot_df$point_order_normalized[i], 1, v[i], alpha = 0.2)
    } else {
      hsv_colors_alpha[i] <- grDevices::hsv(plot_df$point_order_normalized[i], 1, 1, alpha = 0.2)
    }
  }

  # Add to dataframe
  plot_df$hsv_colors <- hsv_colors_alpha

  # Create color palette for the colorbar (no alpha)
  n_colors <- 100
  color_palette <- character(n_colors)
  for(i in 1:n_colors) {
    h_val <- (i-1)/(n_colors-1)
    color_palette[i] <- grDevices::hsv(h_val, 1, 0.8)
  }

  # Create the ggplot
  p <- p +

    geom_path(data = plot_df, ggplot2::aes(x = .data[[cols[1]]],
                                           y = .data[[cols[2]]],
                                           group = 1),
                                           color = hsv_colors_alpha,
                                           size = 1) +
    geom_point(data = plot_df, ggplot2::aes(x = .data[[cols[1]]],
                                            y = .data[[cols[2]]]),
                                            color = hsv_colors_alpha,
                                            size = 2) +
    geom_point(data = plot_df, ggplot2::aes(x = .data[[cols[1]]],
                                            y = .data[[cols[2]]],
                                            color = point_order_normalized),
                                            size = 0,
                                            alpha = 0) +
        # Add color scale for the legend
    ggplot2::scale_color_gradientn(
        colors = color_palette,
        name = "Position",
        guide = ggplot2::guide_colorbar(
          title.position = "top",
          barwidth = 1,
          barheight = 10
          )
        ) +

        # Set plot title and axis labels
    ggplot2::labs(title = title, x = cols[1], y = cols[2]) +
        # Center the title
    ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "right"
        )

  # Apply limits if provided
  p <- p + ggplot2::xlim(xlim) +
    ggplot2::ylim(ylim)

  # Print the plot
  print(p)


  return(p)
}
