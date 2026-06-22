#' Sampling pseudo-absences for the training and testing datasets.
#'
#' \code{paSampling} performs a two-step procedure for uniformly sampling pseudo-absences within the environmental space.
#' In the initial step, a kernel-based filter is utilized to determine the subset of the environmental space that will be subsequently sampled. The kernel-based filter calculates the probability function based on the presence observations, enabling the identification of areas within the environmental space that likely exhibit suitable conditions for the species. To achieve this, a probability threshold value is utilized to assign observations to the corresponding portion of the environmental space. These areas, deemed to have suitable environmental conditions, are excluded from the subsequent uniform sampling process conducted in the second step using the \code{uniformSampling} function, which is internally called.
#' The bandwidth of the kernel can be automatically estimated from the presence observations or directly set by the user, providing flexibility in determining the scope and precision of the filter.
#'
#' Being designed with species distribution models in mind, \code{paSampling} allows collectively sampling pseudo-absences for both the training and testing dataset (optional). In both cases, the user must provide a number of observations that will be sampled in each cell of the sampling grid (\code{n.tr}: points for the training dataset; \code{n.ts}: points for the testing dataset).
#' Note that the optimal resolution of the sampling grid can be found using the \code{optimRes} function. Also, note that the number of pseudo-absences eventually sampled in each cell by the internally-called \code{uniformSampling} function depends on the spatial configuration of the observations within the environmental space. Indeed, in most cases some cells of the sampling grid will be empty (i.e., those at the boundary of the environmental space). For this reason, the number of pseudo-absences returned by \code{paSampling} is likely to be lower than the product between the number of cells of the sampling gird and \code{n.tr}(or \code{n.ts}).
#'
#' @param env.rast A RasterStack, RasterBrick or a SpatRaster object comprising the variables describing the environmental space.
#' @param pres A SpatialPointsDataframe, a SpatVector or an sf object including the presence-only observations of the species of interest.
#' @param thres (double) This value identifies the quantile value used to specify the boundary of the kernel density estimate (default \code{thres=0.75} ). Thus, probability values higher than the threshold should indicate portions of the multivariate space likely associated with presence points.
#' @param H The kernel bandwidth (i.e., the width of the kernel density function that defines its shape) excluding the portion of the environmental space associated with environmental conditions likely suitable for the species. It can be either defined by the user or automatically estimated by \code{paSampling} via \code{ks::Hpi}.
#' @param grid.res (integer) resolution of the sampling grid. The resolution can be arbitrarily selected or defined using the \code{optimRes} function.
#' @param n.tr (integer) number of pseudo-absences for the training dataset to sample in each cell of the sampling grid
#' @param n.ts (integer; optional) number of pseudo-absences for the testing dataset to sample in each cell of the sampling grid. sub.ts argument must be TRUE.
#' @param sub.ts (logical) sample the validation pseudo-absences
#' @param prev (double) prevalence value to be specified instead of n.tr and n.ts
#' @param plot_proc (logical) plot progress of the sampling, default FALSE
#' @param verbose (logical) Print verbose
#' @param dimensions (string vector) specify the pc components to analyse. Must have length >= 2; any number of components is supported (the method is dimension-general, only the compute cost grows with dimension).
#' @param precomputed.pca If in an other step the pca has already been calculated it can be but in here to speed up computation
#' @param nn.based.presence.exclusion Boolean. If true, a candidate pseudo-absence is excluded when it lies in a neighbourhood as dense in presences as the presence core (measured by the mean distance to the nearest presence points). This metric criterion replaces the planar convex hull of the original method and works in any dimension. If false, a kernel-density criterion is used instead (limited to at most 6 dimensions).
#' @param data.based.distance.threshold If true uses the dataset to evaluate realistic distances to the repapped points
#' @param n.samples Number of sample points that are returned
#' @param n.candidates (integer; optional) number of raw uniform proposals drawn per batch over the environmental bounding box. Defaults to \code{grid.res^2 * n.tr} (the legacy two-dimensional budget). Because the realized environment fills an exponentially smaller share of its bounding box as the number of \code{dimensions} grows, high-dimensional runs need a larger value (or a larger \code{n.tr}); the function draws additional batches up to a cap when \code{n.samples} is requested and warns if it cannot reach it.
#' @param dim.correction Dimension correction for the support-membership distance threshold, forwarded to \code{\link{optimalDistanceThresholdNn}}. One of \code{"voronoi"} (default), \code{"simplex"}, \code{"none"}, or a positive numeric multiplier. All equal 1 at two dimensions, so two-dimensional results are unchanged.
#' @importFrom stats na.omit quantile complete.cases runif
#' @return An sf object with the coordinates of the pseudo-absences both in the geographical and environmental space.
#' @export
#'
paSamplingNn <- function (env.rast=NULL, pres = NULL, thres = 0.75, H = NULL, grid.res = 10,
                        n.tr = 5, sub.ts = FALSE, n.ts = 5, prev = NULL, plot_proc = FALSE,
                        verbose = FALSE, dimensions = c("PC1", "PC2"),
                        precomputed.pca = NULL,
                        n.samples = NULL,
                        nn.based.presence.exclusion = TRUE,
                        data.based.distance.threshold = TRUE,
                        n.candidates = NULL,
                        dim.correction = "voronoi") {
  # Input validation
  check_raster_input(env.rast, "env.rast")
  check_spatial_points(pres, "pres", allow_null = TRUE)
  if (!is.character(dimensions) || length(dimensions) < 2) {
    stop("'dimensions' must be a character vector with at least 2 elements", call. = FALSE)
  }
  check_in_range(thres, "thres", min_val = 0, max_val = 1)
  if (!is.numeric(grid.res) || length(grid.res) != 1 || grid.res < 1) {
    stop(paste0("'grid.res' must be a positive number, got ", deparse(grid.res)), call. = FALSE)
  }
  if (!is.numeric(n.tr) || length(n.tr) != 1 || n.tr < 1) {
    stop(paste0("'n.tr' must be a positive number, got ", deparse(n.tr)), call. = FALSE)
  }
  if (!is.null(n.samples) && (!is.numeric(n.samples) || length(n.samples) != 1 || n.samples < 1)) {
    stop(paste0("'n.samples' must be NULL or a positive number, got ", deparse(n.samples)), call. = FALSE)
  }
  if (!is.null(n.candidates) && (!is.numeric(n.candidates) || length(n.candidates) != 1 || n.candidates < 1)) {
    stop(paste0("'n.candidates' must be NULL or a positive number, got ", deparse(n.candidates)), call. = FALSE)
  }
  if (!is.logical(nn.based.presence.exclusion) || length(nn.based.presence.exclusion) != 1) {
    stop(paste0("'nn.based.presence.exclusion' must be a single logical value, got '",
                paste(class(nn.based.presence.exclusion), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.logical(data.based.distance.threshold) || length(data.based.distance.threshold) != 1) {
    stop(paste0("'data.based.distance.threshold' must be a single logical value, got '",
                paste(class(data.based.distance.threshold), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.null(precomputed.pca)) {
    if (!is.list(precomputed.pca) || is.null(precomputed.pca$PCs)) {
      stop("'precomputed.pca' must be a list with a '$PCs' element (result of rastPCA), or NULL", call. = FALSE)
    }
  }

  if (is.null(pres)) {
    message("No species has been supplied. This leads to standard uniform sampling.")
  } else {

    if (is.null(grid.res)) {
      stop("A grid resolution must be provided as length-one integer")
    }
    if (is.null(prev)) {
      estPrev <- round(nrow(pres)/(n.tr * (grid.res^2)), 2)
      message(paste("Estimated prevalence of", estPrev))
    } else {
      n.tr <- (nrow(pres)/prev)/(grid.res^2)
      n.ts <- (nrow(pres)/prev)/(grid.res^2)
      message(paste("User-defined prevalence of", prev))
    }
    if (inherits(env.rast, "BasicRaster")) {
      env.rast <- terra::rast(env.rast)
    }
    if (!inherits(pres, "SpatVector")) {
      occ.vec <- terra::vect(pres)
    } else {
      occ.vec <- pres
    }

  }

  if (!is.null(precomputed.pca)){
    rpc <- precomputed.pca
  } else {
    rpc <- rastPCA(env.rast, stand = TRUE)
  }

  env.data <- env.rast %>%
    as.data.frame(xy = TRUE)

  env.with.pc <- rpc$PCs %>%
    as.data.frame(xy = TRUE) %>%
    na.omit() %>%
    cbind(env.data)

  check_columns_exist(env.with.pc, dimensions, "env.with.pc", "dimensions")

  d <- length(dimensions)
  pc.mat <- as.matrix(env.with.pc[dimensions])

  # --- Uniform proposals over the environmental support (dimension-general) ------
  # The original method built a regular sf grid (sf::st_make_grid) over PC1/PC2 and
  # added per-axis x/y noise. Both are planar-only. We instead draw points directly
  # and uniformly over the per-axis bounding box of the cloud, which is exactly the
  # same proposal distribution in 2D but valid for any number of `dimensions`.
  #
  # Remapping these box-uniform proposals to their nearest real point yields
  # pseudo-absences that are uniform over the SUPPORT in every dimension: a proposal
  # selects a data point iff it lands in that point's Voronoi cell, and the expected
  # Voronoi-cell volume equals 1 / local-density in all d, so the "denser regions
  # hold more points" and "denser regions have smaller cells" effects cancel. The
  # remapping therefore needs NO dimension correction. What does degrade with d is
  # the acceptance rate (below): the cloud fills an exponentially smaller share of
  # its bounding box, so high-dimensional runs are compute-bound by design.
  axis.min  <- apply(pc.mat, 2, min)
  axis.span <- apply(pc.mat, 2, function(x) diff(range(x)))
  # notional per-axis cell spacing, kept so grid.res / n.tr retain their meaning
  step <- axis.span / max(1, grid.res)

  if (is.null(n.candidates)) {
    n.candidates <- max(as.integer(round(grid.res^2 * n.tr)), 1L)
  }
  n.candidates <- as.integer(round(n.candidates))

  # Support-membership distance threshold (dimension-corrected; see
  # optimalDistanceThresholdNn). A proposal is kept only if its nearest real point
  # is within this distance, which rejects proposals that landed in genuine holes
  # outside the support. The sqrt(d/2)-type correction is a no-op at d = 2.
  if (data.based.distance.threshold) {
    distance.threshold <- optimalDistanceThresholdNn(env.data = env.with.pc,
                                                     dimensions = dimensions,
                                                     dim.correction = dim.correction)
  } else {
    distance.threshold <- (max(step) / 2) * .nnThresholdDimFactor(d, dim.correction)
  }

  # Draw, remap and threshold in batches until we have enough unique survivors (or
  # the draw budget is spent). A single batch reproduces the legacy single-pass
  # behaviour when n.samples is not requested.
  target      <- if (!is.null(n.samples)) n.samples else NA_real_
  max.batches <- if (!is.null(n.samples)) 25L else 1L
  kept.index  <- integer(0)
  kept.dist   <- numeric(0)
  n.accepted  <- 0L
  for (b in seq_len(max.batches)) {
    u <- matrix(stats::runif(n.candidates * d), ncol = d)
    u <- sweep(sweep(u, 2, axis.span, `*`), 2, axis.min, `+`)  # scale to [min, max] per axis
    nn <- FNN::get.knnx(pc.mat, u, k = 1)
    keep <- nn$nn.dist[, 1] < distance.threshold
    if (any(keep)) {
      n.accepted <- n.accepted + sum(keep)
      kept.index <- c(kept.index, nn$nn.index[keep, 1])
      kept.dist  <- c(kept.dist,  nn$nn.dist[keep, 1])
      first <- !duplicated(kept.index)           # one row per distinct real point
      kept.index <- kept.index[first]
      kept.dist  <- kept.dist[first]
    }
    if (!is.na(target) && length(kept.index) >= target) break
  }

  if (length(kept.index) == 0L) {
    warning("No proposals fell within the support distance threshold; returning an empty sample. ",
            "In high dimensions the acceptance rate collapses (the realized environment fills an ",
            "exponentially small fraction of its bounding box) - increase 'n.candidates' or 'n.tr'.",
            call. = FALSE)
    return(sf::st_as_sf(env.with.pc[0, , drop = FALSE], coords = c("x", "y")))
  }
  if (!is.na(target) && length(kept.index) < target) {
    message(sprintf(paste0("Only %d unique pseudo-absences survived the support filter after %d ",
                           "batches (requested %d). The acceptance rate falls steeply with ",
                           "dimension; raise 'n.candidates' or 'n.tr' to draw more."),
                    length(kept.index), max.batches, as.integer(target)))
  }

  mapped.sampled.points.filtered <- env.with.pc[kept.index, , drop = FALSE]
  mapped.sampled.points.filtered$distance <- kept.dist
  n.duplicate.maps <- n.accepted - length(kept.index)

  # --- Remove pseudo-absences that fall in the species' presence region ---------
  # The original method excluded points inside the convex hull of the dense presence
  # sub-cloud (sf::st_convex_hull / st_within). A convex hull is planar-only AND, in
  # high dimensions, statistically wrong: it is pinned by a few extreme presences yet
  # encloses vast empty "corner" volume, so it would delete legitimate pseudo-
  # absences far from any presence. We instead apply a local metric criterion in PC
  # space, which reduces to "inside the dense presence blob" in 2D and is valid in
  # any dimension.
  if (!is.null(pres)) {
    virtual.presence.points.pc <- terra::extract(rpc$PCs, occ.vec, bind = TRUE) %>%
      sf::st_as_sf()
    virtual.presence.points.pc <- cbind(sf::st_drop_geometry(virtual.presence.points.pc),
                                        sf::st_coordinates(virtual.presence.points.pc))
    check_columns_exist(virtual.presence.points.pc, dimensions,
                        "virtual.presence.points.pc", "dimensions")
    # Drop presence points that fell on no-data cells; FNN::get.knnx refuses NAs.
    virtual.presence.points.pc <- virtual.presence.points.pc[
      stats::complete.cases(virtual.presence.points.pc[, dimensions, drop = FALSE]), ,
      drop = FALSE]
    if (nrow(virtual.presence.points.pc) < 2) {
      stop("Too few non-NA presence points to compute a presence-density region",
           call. = FALSE)
    }
    pres.mat <- as.matrix(virtual.presence.points.pc[dimensions])
    cand.mat <- as.matrix(mapped.sampled.points.filtered[dimensions])

    if (nn.based.presence.exclusion) {
      if (verbose) message("nearest neighbor based presence exclusion")
      # Presence-density proxy: mean distance to the k nearest OTHER presences.
      k.neighbors <- min(100L, nrow(pres.mat))
      pres.self  <- FNN::get.knnx(pres.mat, pres.mat,
                                  k = min(k.neighbors + 1L, nrow(pres.mat)))
      pres.proxy <- rowMeans(pres.self$nn.dist[, -1, drop = FALSE])   # drop self (col 1)
      # Same statistic at each candidate (distance to its k nearest presences).
      cand.knn   <- FNN::get.knnx(pres.mat, cand.mat, k = k.neighbors)
      cand.proxy <- rowMeans(cand.knn$nn.dist)
      # The presence "core" is the densest (1 - thres) quantile => smallest proxies.
      # Keep candidates that are LESS dense in presences than that core.
      proxy.cut <- quantile(pres.proxy, 1 - thres)
      outside   <- cand.proxy > proxy.cut
    } else {
      if (verbose) message("kernel based presence exclusion")
      if (d > 6L) {
        stop("kernel-based presence exclusion supports at most 6 dimensions (ks::kde); ",
             "use nn.based.presence.exclusion = TRUE for higher-dimensional spaces",
             call. = FALSE)
      }
      if (is.null(H)) H <- ks::Hpi(x = pres.mat)
      pres.kde <- ks::kde(pres.mat, eval.points = pres.mat, H = H)$estimate
      cand.kde <- ks::kde(pres.mat, eval.points = cand.mat, H = H)$estimate
      # 'in' = KDE above the thres quantile of presence-point density; exclude those.
      kde.cut <- quantile(pres.kde, thres)
      outside <- cand.kde <= unname(kde.cut[1])
    }
    sampled.points <- mapped.sampled.points.filtered[outside, , drop = FALSE]
  } else {
    sampled.points <- mapped.sampled.points.filtered
  }

  # Each survivor is already a distinct real point (deduplicated above), so build the
  # output sf directly in geographic space.
  sampled.points.unique <- sf::st_as_sf(sampled.points, coords = c("x", "y"))
  message(paste("\nThere were ", n.duplicate.maps,
                "points that were sampled twice. This indicates undersampling of low density regions or oversampling of the border region.\nThis occures as the probability of beeing close to the same points twice is lower in high denisty regions."))

  if (!is.null(n.samples)){
    sample.indexes <- floor(seq(1, nrow(sampled.points.unique),
                                length.out = min(n.samples, nrow(sampled.points.unique))))
    selected.sampled.points <- sampled.points.unique[sample.indexes, ]
    return(selected.sampled.points)
  }

  return(sampled.points.unique)
}

