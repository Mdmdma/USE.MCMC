#' Function to streamline the generation and sampling of a virtual species
#'
#' @param env.data Raster containing environmental parameters
#' @param n.samples Number of samples
#'
#' @returns List containing the sampled points as a dataframe as well as other nice things that can later be used to plot
#' @export

getVirtualSpeciesPresencePoints <- function(env.data = NULL, n.samples = 0){
  # Create virtual species
  random.sp <- virtualspecies::generateRandomSp(env.data,
                                                convert.to.PA = FALSE,
                                                species.type = "additive",
                                                realistic.sp = TRUE,
                                                plot = FALSE)

  # Reclassify suitability raster using a probability conversion rule
  new.pres <- virtualspecies::convertToPA(x=random.sp,
                                          beta=0.55,
                                          alpha = -0.05, plot = FALSE)

  # Sample true occurrences
  print(n.samples)
  presence.data <- virtualspecies::sampleOccurrences(new.pres,
                                                       n = n.samples, # The number of points to sample
                                                       type = "presence only",
                                                       detection.probability = 1,
                                                       correct.by.suitability = TRUE,
                                                       plot = TRUE)

  # Generate a presence-only data set
  presence.dataset <- presence.data$sample.points[c("x", "y")]
  presence.dataset <- sf::st_as_sf(presence.dataset, coords=c("x", "y"), crs=4326)["geometry"]
  presence.dataset <- terra::vect(presence.dataset)
  presence.data$sample.points <- presence.dataset
  terra::plot(random.sp$suitab.raster, main = "Suitability score given by the model that produced the VS")
  return(presence.data)
}
