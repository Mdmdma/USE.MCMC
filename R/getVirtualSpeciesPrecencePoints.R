#' Function to streamline the generation and sampling of a virtual species
#'
#' @param env.data Raster containing environmental parameters
#' @param n.samples Number of samples
#' @param plot Boolean sets if the function should provide plots of the sampled points and the species suitability.
#' @returns List containing the sampled points as a dataframe as well as other nice things that can later be used to plot
#' @export

getVirtualSpeciesPresencePoints <- function(env.data = NULL, n.samples = 0, plot = FALSE){
  # Input validation
  if (!requireNamespace("virtualspecies", quietly = TRUE)) {
    stop("Package 'virtualspecies' is required for this function but is not installed. ",
         "Install it with install.packages('virtualspecies')", call. = FALSE)
  }
  check_raster_input(env.data, "env.data")
  if (!is.numeric(n.samples) || length(n.samples) != 1 || n.samples < 1 || n.samples != floor(n.samples)) {
    stop(paste0("'n.samples' must be a positive integer, got ", deparse(n.samples)), call. = FALSE)
  }
  if (!is.logical(plot) || length(plot) != 1) {
    stop(paste0("'plot' must be a single logical value, got '",
                paste(class(plot), collapse = "/"), "'"), call. = FALSE)
  }

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
                                                       plot = plot)

  # Generate a presence-only data set
  presence.dataset <- presence.data$sample.points[c("x", "y")]
  presence.dataset <- sf::st_as_sf(presence.dataset, coords=c("x", "y"), crs=4326)["geometry"]
  presence.dataset <- terra::vect(presence.dataset)
  presence.data$sample.points <- presence.dataset
  if (plot){
    terra::plot(random.sp$suitab.raster, main = "Suitability score given by the model that produced the VS")
  }
  return(presence.data)
}
