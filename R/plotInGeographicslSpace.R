plotInGeographicalSpace <- function(presence.map = NULL, presence.points = NULL, absence.points = NULL){
  presence.map <- terra::unwrap(presence.map)
  mfrow = c(1,1)
  plot(presence.map)

  graphics::points(presence.points$geometry, col = "black", pch = 20)
  graphics::points(absence.points$geometry, col = "red", pch = 20)
}
