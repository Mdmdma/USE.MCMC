#' Function that plots the geographical location of points onto a raster
#'
#' @param presence.distribution.raster Raster containing the target distribution
#' @param presence.points Dataframe or list containing a geometry column that contains points to be plotted
#' @param absence.points Dataframe or list containing a geometry column that contains points to be plotted
#' @param minimal Boolean if TRUE removes the labels legend and title
#'
#' @returns NULL, as the function just plots
#' @export
plotInGeographicalSpace <- function(presence.distribution.raster = NULL, presence.points = NULL, absence.points = NULL, minimal = FALSE) {
  # Input validation
  check_raster_input(presence.distribution.raster, "presence.distribution.raster")
  if (is.null(presence.points)) {
    stop("'presence.points' must be provided (got NULL)", call. = FALSE)
  }
  if (!inherits(presence.points, "sf")) {
    stop(paste0("'presence.points' must be an sf object, got '",
                paste(class(presence.points), collapse = "/"), "'"), call. = FALSE)
  }
  if (is.null(absence.points)) {
    stop("'absence.points' must be provided (got NULL)", call. = FALSE)
  }
  if (!inherits(absence.points, "sf")) {
    stop(paste0("'absence.points' must be an sf object, got '",
                paste(class(absence.points), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.logical(minimal) || length(minimal) != 1) {
    stop(paste0("'minimal' must be a single logical value, got '",
                paste(class(minimal), collapse = "/"), "'"), call. = FALSE)
  }
  if (sf::st_crs(presence.points) != sf::st_crs(absence.points)){
    stop("Both points need to be in the same coordinate system", call. = FALSE)
  }

  rs <- terra::classify(presence.distribution.raster, c(-1,0,1))

  # Initialize empty plot
  plot <- ggplot2::ggplot()

  # Add raster layer using tidyterra
  plot<-  plot + tidyterra::geom_spatraster(data = rs)

  plot <- plot + ggplot2::scale_fill_manual(
    values = c("yellow", "purple"), # Assign colors: 0=yellow, 1=violet, 2 (originally NA)=white
    labels = c("Absence", "Presence"),  # Labels for the legend
    name = "VS model",  # Optional: Title for the legend
    na.value = "white",        # Color for NA values in the plot
    na.translate = FALSE,
  )

  # Add presence points with aesthetic mapping for legend
  plot <- plot + ggplot2::geom_sf(
    data = presence.points,
    ggplot2::aes(color = "Presence"),
    size = 1
  )

  # Add absence points with aesthetic mapping for legend
  plot <- plot + ggplot2::geom_sf(
    data = absence.points,
    ggplot2::aes(color = "Pseudo absence"),
    size = 1
  )

  # Add title
  plot <- plot + ggplot2::labs(title = "Geographical position of the points")

  # Add minimal theme
  plot <- plot + ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # Add manual scale for the point legend
  plot <- plot + ggplot2::scale_color_manual(
    name = "Points",
    values = c("Presence" = "black", "Pseudo absence" = "red")
  )
  if (minimal){
    plot <- plot + ggplot2::theme(legend.position = "none",   # remove legend
                         axis.title = ggplot2::element_blank(),   # remove axis titles
                         axis.text = ggplot2::element_blank(),    # remove axis tick labels
                         axis.ticks = ggplot2::element_blank()
                         ) +
      ggplot2::labs(title = NULL, x = NULL, y = NULL)
  }

  return(plot)
}


