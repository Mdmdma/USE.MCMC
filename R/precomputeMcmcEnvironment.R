#' precomputeMcmcEnvironment
#'
#' Runs the *environment-only* portion of [paSamplingMcmc()] once and returns a
#' reusable, serialisable bundle. Everything this function computes depends only
#' on the environmental raster (plus `dimensions`, `seed.number` and the
#' subsampling), and is therefore **constant across presence sets**, so when the
#' same environment is sampled for many species / presence draws, this work can
#' be done a single time and handed back to [paSamplingMcmc()] via its
#' `precomputed.env` argument.
#'
#' Cached steps: the PCA (`rastPCA`), construction of the PC-augmented background
#' data frame, the Gaussian-mixture fit to the environment and its density
#' evaluation over the full background, the proposal covariance matrix, and the
#' nearest-neighbour distance threshold. The species GMM, the combined density
#' function and the MCMC chains themselves remain per-call in
#' [paSamplingMcmc()].
#'
#' The bundle is safe to `saveRDS()` / `readRDS()`: the PC SpatRaster is stored
#' wrapped via [terra::wrap()] (terra objects are external pointers and would not
#' otherwise survive serialisation) and unwrapped on use. This makes the bundle
#' usable across separate R sessions / batch jobs.
#'
#' Reproducibility: the function seeds with `seed.number` and captures the RNG
#' state immediately after the environment block. [paSamplingMcmc()] restores
#' that state before drawing the species model and chains, so results are
#' **identical** whether the environment was computed inline or loaded from a
#' cached bundle.
#'
#' @param env.rast Terra raster containing the environment.
#' @param dimensions Character vector of the PC dimensions to retain (e.g.
#'   `c("PC1", "PC2")`). Must match the `dimensions` later passed to
#'   [paSamplingMcmc()].
#' @param seed.number Seed used for the environment subsampling so the bundle is
#'   reproducible. Must match the `seed.number` later passed to
#'   [paSamplingMcmc()].
#' @param precomputed.pca Optional result of [rastPCA()] (a list with a `$PCs`
#'   element) to reuse instead of recomputing the PCA.
#' @param plot_proc If TRUE, `mclust` draws its fit diagnostics.
#' @param verbose If TRUE, prints progress.
#' @param ... reserved; the deprecated \code{env.data.raster} name (use \code{env.rast}) or any unknown argument raises an informative error.
#'
#' @returns A list (S3 class `mcmc_environment`) with elements `pcs_packed`,
#'   `env.with.pc.sf`, `env.data.cleaned`, `environmental.data.model`,
#'   `environmental.densities`, `covariance.matrix`, `distance.threshold`,
#'   `dimensions`, `seed.number` and `rng_state`. Pass it to
#'   [paSamplingMcmc()]'s `precomputed.env` argument.
#' @export
precomputeMcmcEnvironment <- function(env.rast = NULL,
                                      dimensions = c("PC1", "PC2"),
                                      seed.number = 42,
                                      precomputed.pca = NULL,
                                      plot_proc = FALSE,
                                      verbose = FALSE, ...) {
  # Guide the deprecated raster name and reject unknown arguments (mirrors the
  # samplers' behaviour rather than silently swallowing them via `...`).
  .dots <- list(...)
  if (length(.dots)) {
    if ("env.data.raster" %in% names(.dots)) {
      stop_config("'env.data.raster' has been renamed to 'env.rast'. Pass env.rast = ... instead.")
    }
    stop_config("unknown argument(s) passed to precomputeMcmcEnvironment(): ",
                paste(names(.dots), collapse = ", "))
  }
  # Input validation
  check_raster_input(env.rast, "env.rast")
  if (!is.character(dimensions) || length(dimensions) < 2) {
    stop("'dimensions' must be a character vector with at least 2 elements", call. = FALSE)
  }
  if (!is.numeric(seed.number) || length(seed.number) != 1) {
    stop(paste0("'seed.number' must be a single numeric value, got ", deparse(seed.number)), call. = FALSE)
  }
  if (!is.null(precomputed.pca)) {
    if (!is.list(precomputed.pca) || is.null(precomputed.pca$PCs)) {
      stop("'precomputed.pca' must be a list with a '$PCs' element (result of rastPCA), or NULL", call. = FALSE)
    }
  }

  if (inherits(env.rast, "BasicRaster")) {
    env.rast <- terra::rast(env.rast)
  }

  set.seed(seed.number)

  # Generate the environmental space using PCA
  if (is.null(precomputed.pca)) {
    rpc <- rastPCA(env.rast, stand = TRUE)
  } else {
    rpc <- precomputed.pca
  }

  # Combine environment and PCA layers on the shared raster grid.
  env.with.pc.sf <- terra::as.data.frame(c(env.rast, rpc$PCs), xy = TRUE) %>%
    na.omit() %>%
    sf::st_as_sf(coords = c("x", "y"))

  # subsample env space to speed up the fit
  env.with.pc.sf.subsampled <- env.with.pc.sf[stats::runif(min(nrow(env.with.pc.sf), 2000), 1, nrow(env.with.pc.sf)), ]

  # clean data
  env.data.cleaned <- sf::st_drop_geometry(env.with.pc.sf[dimensions])
  env.data.cleaned.subsampled <- sf::st_drop_geometry(env.with.pc.sf.subsampled[dimensions])

  # environment model
  if (verbose) cat("Fit environmental model \n")
  environmental.data.model <- mclust::densityMclust(env.data.cleaned.subsampled,
                                                    plot = plot_proc,
                                                    verbose = verbose)
  environmental.densities <- stats::predict(environmental.data.model, env.data.cleaned)

  # proposal covariance (env-only)
  covariance.matrix <- stats::cov(sf::st_drop_geometry(env.with.pc.sf)[dimensions])

  # nearest-neighbour distance threshold (env-only)
  distance.threshold <- optimalDistanceThresholdNn(env.data = env.with.pc.sf,
                                                   dimensions = dimensions)

  # Capture the RNG state *after* the environment block so paSamplingMcmc() can
  # restore it and reproduce the exact same species/MCMC stream whether the
  # environment was computed inline or loaded from this cached bundle.
  rng_state <- get(".Random.seed", envir = .GlobalEnv)

  structure(
    list(
      pcs_packed               = terra::wrap(rpc$PCs),
      env.with.pc.sf           = env.with.pc.sf,
      env.data.cleaned         = env.data.cleaned,
      environmental.data.model = environmental.data.model,
      environmental.densities  = environmental.densities,
      covariance.matrix        = covariance.matrix,
      distance.threshold       = distance.threshold,
      dimensions               = dimensions,
      seed.number              = seed.number,
      rng_state                = rng_state
    ),
    class = "mcmc_environment"
  )
}
