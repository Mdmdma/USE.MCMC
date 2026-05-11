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
#' @param burnIn Integer, sets the number of steps per adaptive burnin cycle. If 0 the burnin is skipped
#' @param covariance.correction Integer, sets the inital value of the covariance correction
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
                          burnIn = 1000,
                          covariance.correction = 1,
                          precomputed.pca = NULL,
                          seed.number = 42,
                          n.neighbors.for.statistics = 2, low.end.of.inclueded.points = 100, high.end.of.included.points = 5,
                          environmental.cutof.percentile = 0.001,
                          species.cutoff.threshold = 0.95,
                          plot_proc = FALSE,
                          num.chains = 1, num.cores = 1) {
  # Input validation
  check_raster_input(env.data.raster, "env.data.raster")
  check_spatial_points(pres, "pres")
  if (!is.numeric(n.samples) || length(n.samples) != 1 || n.samples < 1) {
    stop(paste0("'n.samples' must be a positive number, got ", deparse(n.samples)), call. = FALSE)
  }
  if (!is.numeric(chain.length) || length(chain.length) != 1 || chain.length < 1) {
    stop(paste0("'chain.length' must be a positive number, got ", deparse(chain.length)), call. = FALSE)
  }
  if (!is.character(dimensions) || length(dimensions) < 2) {
    stop("'dimensions' must be a character vector with at least 2 elements", call. = FALSE)
  }
  if (!is.numeric(burnIn) || length(burnIn) != 1 || burnIn < 0) {
    stop(paste0("'burnIn' must be a non-negative number, got ", deparse(burnIn)), call. = FALSE)
  }
  if (!is.numeric(covariance.correction) || length(covariance.correction) != 1 || covariance.correction <= 0) {
    stop(paste0("'covariance.correction' must be a positive number, got ", deparse(covariance.correction)), call. = FALSE)
  }
  if (!is.null(precomputed.pca)) {
    if (!is.list(precomputed.pca) || is.null(precomputed.pca$PCs)) {
      stop("'precomputed.pca' must be a list with a '$PCs' element (result of rastPCA), or NULL", call. = FALSE)
    }
  }
  if (!is.numeric(seed.number) || length(seed.number) != 1) {
    stop(paste0("'seed.number' must be a single numeric value, got ", deparse(seed.number)), call. = FALSE)
  }
  if (!is.numeric(num.chains) || length(num.chains) != 1 || num.chains < 1) {
    stop(paste0("'num.chains' must be a positive integer, got ", deparse(num.chains)), call. = FALSE)
  }
  if (!is.numeric(num.cores) || length(num.cores) != 1 || num.cores < 1) {
    stop(paste0("'num.cores' must be a positive integer, got ", deparse(num.cores)), call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop(paste0("'verbose' must be a single logical value, got '", paste(class(verbose), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.logical(plot_proc) || length(plot_proc) != 1) {
    stop(paste0("'plot_proc' must be a single logical value, got '", paste(class(plot_proc), collapse = "/"), "'"), call. = FALSE)
  }
  check_in_range(environmental.cutof.percentile, "environmental.cutof.percentile", min_val = 0, max_val = 1)
  check_in_range(species.cutoff.threshold, "species.cutoff.threshold", min_val = 0, max_val = 1)

  if (inherits(env.data.raster, "BasicRaster")) {
    env.data.raster <- terra::rast(env.data.raster)
  }
  set.seed(seed.number)

  # Generate the environmental space using PCA
  if (is.null(precomputed.pca)){
    rpc <- rastPCA(env.data.raster, stand = TRUE)
  } else {
    rpc <- precomputed.pca
  }

  # Combine environment and PCA layers on the shared raster grid; avoids the
  # expensive sf::st_join over tens-of-thousands of point geometries that the
  # earlier pipeline performed.
  env.with.pc.sf <- terra::as.data.frame(c(env.data.raster, rpc$PCs), xy = TRUE) %>%
    na.omit() %>%
    sf::st_as_sf(coords = c("x", "y"))

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
  species.cutoff.threshold <- stats::quantile(species.densities, species.cutoff.threshold)

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
  results.computation <- parallel::mclapply(1:num.chains, function(interator) {
  # sample points
    sampled.points <- mcmcSampling(dataset = env.with.pc.sf,
                                   dimensions = dimensions,
                                   n.sample.points = chain.length,
                                   proposalFunction = proposalFunction,
                                   densityFunction = densityFunction,
                                   burnIn = burnIn,
                                   covariance.correction = covariance.correction,
                                   verbose = verbose)
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


