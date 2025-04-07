getVirtualSpeciesPresencePoints <- function(environemtalData = NULL, n.samples = 0){
  # Create virtual species
  random.sp <- virtualspecies::generateRandomSp(environemtalData,
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
  presence.points <- virtualspecies::sampleOccurrences(new.pres,
                                                       n = 100, # The number of points to sample
                                                       type = "presence-absence",
                                                       sample.prevalence = 0.99,
                                                       detection.probability = 1,
                                                       correct.by.suitability = TRUE,
                                                       plot = TRUE)

  # Generate a presence-only data set
  myPres <- presence.points$sample.points[which(presence.points$sample.points$Observed==1), c("x", "y", "Observed")]
  myPres <- sf::st_as_sf(myPres, coords=c("x", "y"), crs=4326)["geometry"]
  myPres <- terra::vect(myPres)

}
