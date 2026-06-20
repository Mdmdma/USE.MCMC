#' paSamplingMcmc
#' is a near drop in replacement for paSampling from the original USE package, that allows to perform a Gaussian mixture based pseudo
#' absence sampling using a markov. In a first step a density function is constructed using a GMM fitted to the environment as a
#' limit to the sampling space and a GMM fitted on the target species as a way to evade regions associated with the presence.
#'
#' @param env.data.raster Terra raster containing the environment. Required unless a `precomputed.env` bundle is supplied, in which case it is optional (the environment is already built and presence PC scores are read from the cached PC rasters).
#' @param pres Sf dataframe containing the presence locations
#' @param n.samples number of samples that should be put out
#' @param chain.length number of points that are sampled for the chain
#' @param verbose If true the function gives updates on the current state of the chain
#' @param dimensions vector containg the names of the dimensions that should be included
#' @param burnIn Integer, number of Robbins-Monro burn-in adaptation steps performed before sampling. During each step the proposal scale is adjusted toward target acceptance 0.234 (Roberts/Rosenthal 2009). Set to 0 to skip adaptation and start sampling immediately at the user-supplied `covariance.correction`.
#' @param covariance.correction Integer, sets the inital value of the covariance correction
#' @param precomputed.pca If rastPCA has already been evoked, it the result of it can be passed here to not recompute
#' @param precomputed.env Optional environment bundle returned by [precomputeMcmcEnvironment()]. The environment fit and PCA are constant across presence sets, so when sampling many species on the same environment they can be computed once and passed here to skip the PCA, the environmental GMM fit/density evaluation, the proposal covariance and the distance-threshold computation. The bundle is `saveRDS()`-safe and can be reused across R sessions / batch jobs. When supplied it takes precedence over `precomputed.pca`; its captured RNG state is restored so results match the inline (uncached) path exactly.
#' @param seed.number seednumber used to get repeatable results
#' @param n.neighbors.for.statistics number of neighbors used to calculate the maximal sensible distance to real points that should be included
#' @param low.end.of.inclueded.points Sets the range of points included in the threshold computation
#' @param high.end.of.included.points Sets the range of points included in the threshold computation
#' @param environmental.cutof.percentile sets the percentile of the environment GMM that is excluded from the space that can be visited by the chain
#' @param species.cutoff.threshold sets the percentile of the species presence GMM that is included in the space that can be visited by the chain
#' @param plot_proc If true the function returns plots the progress
#' @param num.chains Number of chains from which samples should be picked
#' @param num.cores Number of cores available for parallelization of the multi-chain computation
#' @param engine One of `"auto"` (default), `"R"`, or `"cpp"`. `"auto"` picks the C++ inner loop when both the internal density and proposal functions are built by `mclustDensityFunction()` and `addHighDimGaussian()` (which is the case here) and falls back to the R loop otherwise. `"cpp"` forces the C++ path. `"R"` forces the pure-R reference loop.
#'
#' @returns dataframe containing the sampled points
#' @export

paSamplingMcmc <- function (env.data.raster=NULL, pres = NULL, n.samples = 300, chain.length = 10000,
                          verbose = FALSE, dimensions = c("PC1", "PC2"),
                          burnIn = 1000,
                          covariance.correction = 1,
                          precomputed.pca = NULL,
                          precomputed.env = NULL,
                          seed.number = 42,
                          n.neighbors.for.statistics = 2, low.end.of.inclueded.points = 100, high.end.of.included.points = 5,
                          environmental.cutof.percentile = 0.001,
                          species.cutoff.threshold = 0.95,
                          plot_proc = FALSE,
                          num.chains = 1, num.cores = 1,
                          engine = c("auto", "R", "cpp")) {
  engine <- match.arg(engine)
  # Input validation. env.data.raster is only needed to build the environment
  # (PCA + GMM fit); when a precomputed.env bundle is supplied that work is
  # already done and the raster is optional.
  if (is.null(precomputed.env) || !is.null(env.data.raster)) {
    check_raster_input(env.data.raster, "env.data.raster")
  }
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
  if (!is.null(precomputed.env)) {
    # rng_state is required: without it the restore below is skipped and the
    # species/MCMC stream would silently diverge from the uncached path.
    required.env.fields <- c("pcs_packed", "env.with.pc.sf", "env.data.cleaned",
                             "environmental.data.model", "environmental.densities",
                             "covariance.matrix", "distance.threshold", "dimensions",
                             "rng_state")
    if (!is.list(precomputed.env) || !all(required.env.fields %in% names(precomputed.env))) {
      stop_config("'precomputed.env' must be the list returned by precomputeMcmcEnvironment() (or NULL). Missing fields: ",
                  paste(setdiff(required.env.fields, names(precomputed.env)), collapse = ", "))
    }
    if (!identical(precomputed.env$dimensions, dimensions)) {
      stop_config("'precomputed.env' was built for dimensions (",
                  paste(precomputed.env$dimensions, collapse = ", "),
                  ") but 'dimensions' is (", paste(dimensions, collapse = ", "),
                  "). Recompute the bundle for the requested dimensions.")
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
  # Seeds the inline environment block below; overridden by the RNG-state restore
  # when a precomputed.env bundle is supplied (kept as a fallback otherwise).
  set.seed(seed.number)

  # Environment block (PCA + environmental GMM fit/density + proposal covariance
  # + distance threshold). This is constant across presence sets, so it is
  # factored into precomputeMcmcEnvironment(): computed inline here when no cache
  # is supplied, or reused verbatim from `precomputed.env`. Sharing the one code
  # path guarantees the cached and uncached results are identical.
  if (is.null(precomputed.env)) {
    precomputed.env <- precomputeMcmcEnvironment(env.data.raster = env.data.raster,
                                                 dimensions = dimensions,
                                                 seed.number = seed.number,
                                                 precomputed.pca = precomputed.pca,
                                                 plot_proc = plot_proc,
                                                 verbose = verbose)
  }

  # Restore the RNG state captured right after the environment block so the
  # species model and MCMC chains draw the same random stream whether the
  # environment was computed inline (above) or loaded from a cached bundle. For
  # the inline path this is a no-op (state is already there).
  if (!is.null(precomputed.env$rng_state)) {
    assign(".Random.seed", precomputed.env$rng_state, envir = .GlobalEnv)
  }

  rpc <- list(PCs = terra::unwrap(precomputed.env$pcs_packed))
  env.with.pc.sf <- precomputed.env$env.with.pc.sf
  env.data.cleaned <- precomputed.env$env.data.cleaned
  environmental.data.model <- precomputed.env$environmental.data.model
  environmental.densities <- precomputed.env$environmental.densities
  covariance.matrix <- precomputed.env$covariance.matrix
  distance.threshold <- precomputed.env$distance.threshold
  # Threshold stays a live tunable: cheap quantile over the cached densities.
  environment.threshold <- stats::quantile(environmental.densities, environmental.cutof.percentile)

  # sample species model

  virtual.presence.points <- pres

  # Presence PC scores come from the (cached or freshly computed) PC rasters
  # alone; the original env layers are not needed here, so env.data.raster is
  # optional once a precomputed.env bundle is supplied.
  virtual.presence.points.pc <- terra::extract(rpc$PCs, virtual.presence.points, bind = TRUE) %>%
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


  # # set sampling parameters (covariance.matrix comes from precomputed.env)
  covariance.scaling <-0.075
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
                                   verbose = verbose,
                                   engine = engine)
  }, mc.cores = min(num.chains, num.cores))

  sampled.points <- do.call(rbind, results.computation)

  mapped.sampled.point.locations <- FNN::get.knnx(env.data.cleaned[dimensions], sampled.points[dimensions],k = 1)
  mapped.sampled.points <- env.with.pc.sf[mapped.sampled.point.locations$nn.index,]
  mapped.sampled.points$density <- sampled.points$density
  mapped.sampled.points$distance <- mapped.sampled.point.locations$nn.dist

  # distance.threshold comes from precomputed.env (env-only quantity)
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


