#' paSamplingMcmc
#' is a near drop in replacement for paSampling from the original USE package, that allows to perform a Gaussian mixture based pseudo
#' absence sampling using a markov. In a first step a density function is constructed using a GMM fitted to the environment as a
#' limit to the sampling space and a GMM fitted on the target species as a way to evade regions associated with the presence.
#'
#' @param env.data.raster Terra raster containing the environment
#' @param pres Sf dataframe containing the presence locations
#' @param n.samples number of samples that should be put out
#' @param chain.length number of points that are sampled for the chain
#' @param verbose If true the function gives updates on the current state of the chain
#' @param dimensions vector containg the names of the dimensions that should be included
#' @param burn.in If False the burnin is skipped
#' @param precomputed.pca If rastPCA has already been evoked, it the result of it can be passed here to not recompute
#' @param seed.number seednumber used to get repeatable results
#' @param n.neighbors.for.statistics number of neighbors used to calculate the maximal sensible distance to real points that should be included
#' @param low.end.of.inclueded.points Sets the range of points included in the threshold computation
#' @param high.end.of.included.points Sets the range of points included in the threshold computation
#' @param environmental.cutof.percentile sets the percentile of the environment GMM that is excluded from the space that can be visited by the chain
#' @param species.cutoff.threshold sets the percentile of the species presence GMM that is included in the space that can be visited by the chain
#' @param plot_proc If true the function returns plots the progress
#' @param num.chains Number of chains from which samples should be picked
#' @param num.cores Number of cores available for parallelization of the multi-chain computation
#'
#' @returns dataframe containing the sampled points
#' @export

paSamplingMcmc <- function (env.data.raster=NULL, pres = NULL, n.samples = 300, chain.length = 10000,
                          verbose = FALSE, dimensions = c("PC1", "PC2"),
                          burn.in = TRUE,
                          precomputed.pca = NULL,
                          seed.number = 42,
                          n.neighbors.for.statistics = 2, low.end.of.inclueded.points = 100, high.end.of.included.points = 5,
                          environmental.cutof.percentile = 0.001,
                          species.cutoff.threshold = 0.95,
                          plot_proc = FALSE,
                          num.chains = 1, num.cores = 1) {

  env.data.sf <- env.data.raster %>%
    as.data.frame(xy = TRUE) %>%
    sf::st_as_sf(coords = c("x", "y"))

  # fixing the
  set.seed(seed.number)

  # Generate the environmental space using PCA
  if (is.null(precomputed.pca)){
    rpc <- rastPCA(env.data.raster,  stand = TRUE)
  } else {
    rpc <- precomputed.pca
  }

  # Attaching the data in the PCA coordinates
  env.with.pc.sf <- rpc$PCs %>%
    as.data.frame(xy = TRUE) %>%
    na.omit() %>%
    sf::st_as_sf(coords = c("x", "y")) %>%
    sf::st_join(env.data.sf)

  # subsample env space to speed up the process
  env.with.pc.sf.subsampled <- env.with.pc.sf[stats::runif(min(nrow(env.with.pc.sf), 2000) , 1, nrow(env.with.pc.sf)),]

  # clean data
  env.data.cleaned <- sf::st_drop_geometry(env.with.pc.sf[dimensions])
  env.data.cleaned.subsampled <- sf::st_drop_geometry(env.with.pc.sf.subsampled[dimensions])

  # environment model
  if (verbose) cat("Fit environmental model \n")
  environmental.data.model <- mclust::densityMclust(env.data.cleaned.subsampled,
                                                    plot = plot_proc,
                                                    verbose = verbose )
  summary(environmental.data.model)
  environmental.densities <- mclust::predict.densityMclust(environmental.data.model, env.data.cleaned)
  environment.threshold <- stats::quantile(environmental.densities, environmental.cutof.percentile)

  # sample species model

  virtual.presence.points <- pres

  env.data.raster.with.pc <- c(env.data.raster, rpc$PCs)
  virtual.presence.points.pc <- terra::extract(env.data.raster.with.pc, virtual.presence.points, bind = TRUE) %>%
    sf::st_as_sf()
  if (verbose) cat("Fit presence model \n")
  species.model = mclust::densityMclust(sf::st_drop_geometry(virtual.presence.points.pc[dimensions]),
                                        plot = plot_proc,
                                        verbose = verbose)
  summary(species.model)
  species.densities <- species.model$density
  species.cutoff.threshold <- stats::quantile(species.densities, 0.9)

  #density Function
  densityFunction <- mclustDensityFunction(env.model = environmental.data.model,
                                           species.model = species.model,
                                           dim = dimensions,
                                           threshold = environment.threshold,
                                           species.cutoff.threshold = species.cutoff.threshold)


  # # set sampling parameters
  covariance.scaling <-0.075
  covariance.matrix <- stats::cov(sf::st_drop_geometry(env.with.pc.sf)[dimensions])
  proposalFunction <- addHighDimGaussian(cov.mat =covariance.scaling * covariance.matrix,
                                         dim = length(dimensions))

  # Set up for multiple chains that sampled from
  results.computation <- list()
  results.computation <- mclapply(1:num.chains, function(interator) {
  # sample points
    sampled.points <- mcmcSampling(dataset = env.with.pc.sf,
                                   dimensions = dimensions,
                                   n.sample.points = chain.length,
                                   proposalFunction = proposalFunction,
                                   densityFunction = densityFunction,
                                   burnIn = burn.in,
                                   verbose = TRUE)
  }, mc.cores = min(num.chains, num.cores))

  sampled.points <- do.call(rbind, results.computation)

  mapped.sampled.point.locations <- FNN::get.knnx(env.data.cleaned[dimensions], sampled.points[dimensions],k = 1)
  mapped.sampled.points <- env.with.pc.sf[mapped.sampled.point.locations$nn.index,]
  mapped.sampled.points$density <- sampled.points$density
  mapped.sampled.points$distance <- mapped.sampled.point.locations$nn.dist

  distance.threshold <- optimalDistanceThresholdNn(env.data = env.with.pc.sf,
                                                   dimensions = dimensions)
  filtered.mapped.sampled.points <- mapped.sampled.points[mapped.sampled.points$distance < distance.threshold, ]

  sample.indexes <- floor(seq(1, nrow(filtered.mapped.sampled.points), length.out = min(n.samples * 2, nrow(filtered.mapped.sampled.points))))
  filtered.mapped.sampled.points.subselected <- filtered.mapped.sampled.points[sample.indexes, ]

  # The selection is done in two steps to be able to return the exact amount of desired points
  # without duplicates. If we remove duplicatates before selecting points, we would undersample
  # low density regions

  filtered.mapped.sampled.points.subselected.unique <- filtered.mapped.sampled.points.subselected[
    !duplicated(filtered.mapped.sampled.points.subselected[dimensions[1]][]),]
  if (verbose){
    message(paste("\nThere were ", nrow(filtered.mapped.sampled.points.subselected) - nrow(filtered.mapped.sampled.points.subselected.unique),
                  "points that were sampled twice. A high number indicates undersampling in low density regions or oversampling at the border regions.
                In case of the first reduce the number of samples, in case of the later there is little to be done"))
  }

  sample.indexes <- floor(seq(1, nrow(filtered.mapped.sampled.points.subselected.unique),
                              length.out = min(n.samples, nrow(filtered.mapped.sampled.points.subselected.unique))))
  final.sampled.points <- filtered.mapped.sampled.points.subselected.unique[sample.indexes, ]

  # TODO check where the coordinate system gets lost
  sf::st_crs(final.sampled.points) <- 4326

  return(final.sampled.points)
}


