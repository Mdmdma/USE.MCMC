#' plotDensityLines
#' enables the plotting of density function as well as a trace of a chain made of points in a dataframe
#' on this density surface. In addition the species that was used to generate the density model can be supplied
#' to verify the performance of the model. If that is the case the supplied dataset will be plotted as points
#' instead of as a chain. This function has some issues with the check of devtools
#'
#'
#' @param dataset Dataframe specifying the dataset to be plotted as a line or points respectively
#' @param xlim x-limit of the plot
#' @param ylim y-limit of the plot
#' @param title title of the plot
#' @param lines boolean that tells the function if the lines should be plotted. If false the dataset will be ploted as points
#' @param cols specifies the columns of the dataframe that are used for the density computation. It has to match to the columns used to build the density model.
#' @param density boolean that tells the function if the density should be plotted
#' @param species dataframe containing the simulated species, it has to contain the column names of cols
#' @param densityFunction a function that can take a dataframe with the columns given in cols and retruns the density at that point
#' @param resolution int sets the resolution of the denisity grid
#' @param minimal boolean if true removes titel, label, legend and axis
#'
#' @returns the greated plot
#' @export
#'
# TODO Fix issues with devtools check
plotDensityLines <- function(dataset, xlim = c(0,1), ylim = c(0,1),
                             title = "Connected Data Points",
                             lines = FALSE, cols = NULL,
                             density = FALSE, species = NULL, densityFunction = alwaysOne, resolution = 10,
                             minimal = FALSE) {
  p <- ggplot2::ggplot()

  # plot density
  if (density){
    if(is.null(dataset)) stop("at least one line of data has to be supplied as a template for the density function")

    # Create a grid matrix
    x_seq <- seq(xlim[1], xlim[2], length.out = resolution)
    y_seq <- seq(ylim[1], ylim[2], length.out = resolution)
    grid <- expand.grid(x = x_seq, y = y_seq)
    colnames(grid) <- c("PC1", "PC2")

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

    # Add filled tiles with explicit data and aesthetics
    p <- p + ggplot2::geom_tile(mapping = ggplot2::aes(x = PC1, y = PC2, fill = density), data = grid)
    p <- p + ggplot2::scale_fill_viridis_c(name = "Density")

    # Add contour lines with explicit data and aesthetics
    p <- p + ggplot2::geom_contour(mapping = ggplot2::aes(x = PC1, y = PC2, z = density), data = grid, color = "black", alpha = 0.3, linewidth = 0.5)

    if (!lines){
      # Add points from dataset if available
      if (nrow(dataset) > 1) {
        p <- p + ggplot2::geom_point(data = dataset,
                                     mapping = ggplot2::aes(x = PC1, y = PC2, color = "Samples"),
                                     inherit.aes = FALSE)
      }

      # Add species points in red if provided
      if (!is.null(species)) {
        p <- p + ggplot2::geom_point(data = species,
                                     mapping = ggplot2::aes(x = PC1, y = PC2, color = "Species"),
                                     inherit.aes = FALSE)
      }

      # Add color scale for points
      p <- p + ggplot2::scale_color_manual(name = "Points",
                                           values = c("Samples" = "red", "Species" = "black"))
    }


  }
  #plot lines if wanted
  if (lines) {
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
      }
    }
    if (density) {
      alpha <- 0.6
    }  else alpha <- 0.25
    # Create HSV colors for points and lines - with alpha
    hsv_colors_alpha <- character(nrow(plot_df))
    for(i in 1:nrow(plot_df)) {
      if(is.numeric(v) && length(v) > 1) {
        hsv_colors_alpha[i] <- grDevices::hsv(plot_df$point_order_normalized[i], 1, v[i], alpha = alpha)
      } else {
        hsv_colors_alpha[i] <- grDevices::hsv(plot_df$point_order_normalized[i], 1, 1, alpha = alpha)
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

      ggplot2::geom_path(data = plot_df, ggplot2::aes(x = .data[[cols[1]]],
                                             y = .data[[cols[2]]],
                                             group = 1),
                                             color = hsv_colors_alpha,
                                             size = 1) +
      ggplot2::geom_point(data = plot_df, ggplot2::aes(x = .data[[cols[1]]],
                                              y = .data[[cols[2]]]),
                                              color = hsv_colors_alpha,
                                              size = 2) +
      ggplot2::geom_point(data = plot_df, ggplot2::aes(x = .data[[cols[1]]],
                                              y = .data[[cols[2]]],
                                              color = point_order_normalized),
                                              size = 0,
                                              alpha = 0) +
      # Add color scale for the legend
      ggplot2::scale_color_gradientn(colors = color_palette,
                                     name = "Position",
                                     guide = ggplot2::guide_colorbar(title.position = "top",
                                                                     barwidth = 1,
                                                                     barheight = 6),
                                                                     labels = c(1,
                                                                                nrow(plot_df)* 0.25,
                                                                                nrow(plot_df) *0.5,
                                                                                nrow(plot_df) *0.75,
                                                                                nrow(plot_df))
                                                                     ) +

      # Set plot title and axis labels
      ggplot2::labs(title = title,
                    x = cols[1],
                    y = cols[2]) +
      # Center the title
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                                        legend.position = "right")
  }


   # Set plot limits and labels
  p <- p + ggplot2::xlim(xlim) +
           ggplot2::ylim(ylim) +
           ggplot2::labs(x = "PC1", y = "PC2",
                         title = title,
                         fill = "Density")

  # Add theme elements
  p <- p + ggplot2::theme_minimal() +
           ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # Print the plot
  if (minimal){
    p <- p + ggplot2::theme(legend.position = "none",   # remove legend
                         axis.title = ggplot2::element_blank(),   # remove axis titles
                         axis.text = ggplot2::element_blank(),    # remove axis tick labels
                         axis.ticks = ggplot2::element_blank()
    ) +
      ggplot2::labs(title = NULL, x = NULL, y = NULL)
  }

  return(p)
}
