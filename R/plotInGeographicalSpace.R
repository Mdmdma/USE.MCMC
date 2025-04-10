#' Function that plots the geographical location of points onto a raster
#'
#' @param presence.distribution.raster Raster containing the target distribution
#' @param presence.points Dataframe or list containing a geometry cloumn that contains points to be plotted
#' @param absence.points Dataframe or list containing a geometry cloumn that contains points to be plotted
#'
#' @returns NULL, as the function just plots
#' @export
plotInGeographicalSpace <- function(presence.distribution.raster = NULL, presence.points = NULL, absence.points = NULL){
  # Create plot with no legend
  terra::plot(presence.distribution.raster, main = "Geographical position of the points", legend = FALSE)

  # Add points
  graphics::points(presence.points$geometry, col = "black", pch = 20)
  graphics::points(absence.points$geometry, col = "red", pch = 20)

  # Force legend to appear by setting xpd=TRUE (allows plotting outside figure region)
  graphics::par(xpd = TRUE)
  graphics::legend(
    "topright",
    legend = c("Presence", "Absence"),
    pch = c(20, 20),
    col = c("black", "red"),
    bty = "n",
    inset = c(0.01, 0.01)  # Small inset to avoid edge issues
  )
  graphics::par(xpd = FALSE)  # Reset parameter
}
