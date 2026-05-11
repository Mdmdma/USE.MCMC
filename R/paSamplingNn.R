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
#' @param dimensions (string vector) specify the pc components to analyse. Has to have length 2
#' @param precomputed.pca If in an other step the pca has already been calculated it can be but in here to speed up computation
#' @param nn.based.presence.exclusion Boolean, if true uses nn to select the quantile with closest neighbors and computes the convex hull around them. All points laying inside this convex hull are excludet.
#' @param data.based.distance.threshold If true uses the dataset to evaluate realistic distances to the repapped points
#' @param n.samples Number of sample points that are returned
#' @importFrom stats na.omit quantile
#' @return An sf object with the coordinates of the pseudo-absences both in the geographical and environmental space.
#' @export
#'
paSamplingNn <- function (env.rast=NULL, pres = NULL, thres = 0.75, H = NULL, grid.res = 10,
                        n.tr = 5, sub.ts = FALSE, n.ts = 5, prev = NULL, plot_proc = FALSE,
                        verbose = FALSE, dimensions = c("PC1", "PC2"),
                        precomputed.pca = NULL,
                        n.samples = NULL,
                        nn.based.presence.exclusion = TRUE,
                        data.based.distance.threshold = TRUE) {
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

  env.with.pc.fs <- sf::st_as_sf(env.with.pc, coords = dimensions)

  # compute the grid cells
  grid <- sf::st_make_grid(env.with.pc.fs, n = c(grid.res)) %>%
    sf::st_centroid() %>%
    sf::st_coordinates() %>%
    as.data.frame()
  names(grid) = dimensions
  # repeat the grid so to get the desired number of points per cell
  grid.repeated <- do.call(rbind, replicate(n.tr, grid, simplify = FALSE))

  # add uniform noise to simulate sampling a random location in each cell.
  # Use the grid extent + resolution to derive the per-axis cell spacing; this is
  # robust to st_make_grid's cell ordering (column-major vs row-major).
  step.scaling = 1
  step.x <- (diff(range(grid[, 1])) / max(1, grid.res - 1)) / step.scaling
  step.y <- (diff(range(grid[, 2])) / max(1, grid.res - 1)) / step.scaling

  noise.x <- stats::runif(nrow(grid.repeated), -step.x, step.x)
  noise.y <- stats::runif(nrow(grid.repeated), -step.y, step.y)
  noise <- data.frame(noise.x, noise.y)
  names(noise) <- dimensions
  grid.noisy <- grid.repeated + noise

  # map the sampled points to real points in the environments
  mapped.sampled.point.data <- FNN::get.knnx(env.with.pc[dimensions], grid.noisy, k = 1)
  mapped.sampled.points <- env.with.pc[mapped.sampled.point.data$nn.index,]
  mapped.sampled.points$distance <- mapped.sampled.point.data$nn.dist

  if (data.based.distance.threshold) {
    distance.threshold <- optimalDistanceThresholdNn(env.data = env.with.pc,
                                                     dimensions = dimensions)
  } else {
    distance.threshold <- max(step.y, step.x) / 2
  }
  mapped.sampled.points.filtered <- mapped.sampled.points[mapped.sampled.points$distance < distance.threshold, ]

  # Remove points that are located in the region that is associated with the target species
  # Analog to the paper the convex hull is computed, but the we just remove the points that layed in this area, instead of
  # excluding them from the possible dataset. This is done as it could happen that the maximal grid resolution could change excluding this area.
  if (!is.null(pres)) {
    # use distance to the 100 closest other presence locations as a proxy for density instead of the kernel.
    if (nn.based.presence.exclusion) {
      if (verbose) message("nearest neigbor based presence exclusion")
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
        stop("Too few non-NA presence points to compute a presence-density convex hull",
             call. = FALSE)
      }

      k.neighbors <- min(100, nrow(virtual.presence.points.pc))
      get.nearest.neighbors <- FNN::get.knnx(virtual.presence.points.pc[dimensions],
                                             virtual.presence.points.pc[dimensions],
                                             k = k.neighbors)
      average.distance.to.neighbors <- rowMeans(get.nearest.neighbors$nn.dist)
      idx <- which(average.distance.to.neighbors <= quantile(average.distance.to.neighbors, 1 - thres))
      points.with.closest.neighbors <- virtual.presence.points.pc[idx, ]
      convex.hull <- sf::st_as_sf(points.with.closest.neighbors, coords = dimensions) %>%
        sf::st_union() %>%
        sf::st_convex_hull()
      point.data.sf <- sf::st_as_sf(mapped.sampled.points.filtered, coords = dimensions)
      outside.of.the.region.with.presence <- point.data.sf[!sf::st_within(point.data.sf, convex.hull, sparse = FALSE), ]
      sampled.points <- sf::st_coordinates(outside.of.the.region.with.presence)
      colnames(sampled.points) <- dimensions
      sampled.points <- cbind(sf::st_drop_geometry(outside.of.the.region.with.presence), sampled.points)
    } else {
      if (verbose) message("kernel based presence exclusion")
      id_rast <- terra::rast(vals = 1:terra::ncell(env.rast),
                             names = "myID",
                             extent = terra::ext(env.rast),
                             nrows = terra::nrow(env.rast),
                             ncols = terra::ncol(env.rast),
                             crs = terra::crs(env.rast))
      abio.st <- terra::as.data.frame(c(id_rast, env.rast))
      dt <- terra::as.data.frame(c(id_rast, rpc$PCs[[dimensions]]), xy = TRUE)
      PC12occ <- terra::extract(id_rast, occ.vec, cells = FALSE, df = FALSE, ID = FALSE)[, 1]
      PC12ex <- na.omit(data.frame(dt, PA = ifelse(dt$myID %in% PC12occ, 1, 0)))
      if (is.null(H)) {
        H <- ks::Hpi(x = PC12ex[, dimensions])
      }

      estimate <- data.frame(KDE = ks::kde(PC12ex[PC12ex$PA == 1, dimensions],
                                           eval.points = PC12ex[PC12ex$PA == 1, dimensions],
                                           h = H)$estimate,
                             PC12ex[PC12ex$PA == 1, c(dimensions, "myID", "PA")])

      quantP <- quantile(estimate[, "KDE"], thres)
      estimate$percP <- ifelse(estimate$KDE <= unname(quantP[1]), "out", "in")
      point.data <- merge(x = PC12ex, y = estimate[estimate$PA == 1, c("myID", "percP")],
                          by = "myID", all.x = TRUE)
      point.data$percP <- ifelse(is.na(point.data$percP), "pabs", point.data$percP)
      chull <- sf::st_as_sf(subset(point.data, point.data$percP == "in", select = dimensions),
                            coords = dimensions)
      chull <- sf::st_union(chull)
      chull <- sf::st_convex_hull(chull)
      point.data.sf <- sf::st_as_sf(mapped.sampled.points.filtered, coords = dimensions)
      outside.of.the.region.with.presence <- point.data.sf[!sf::st_within(point.data.sf, chull, sparse = FALSE), ]
      sampled.points <- sf::st_coordinates(outside.of.the.region.with.presence)
      colnames(sampled.points) <- dimensions
      sampled.points <- cbind(sf::st_drop_geometry(outside.of.the.region.with.presence), sampled.points)
    }
  } else {
    sampled.points <- mapped.sampled.points.filtered
  }
  # Reorganize the columns of the sf object so that the output is more convenient to use

  # select only unique points
  sampled.points.unique <- sampled.points[!duplicated(sampled.points[[dimensions[1]]]), ] %>%
    sf::st_as_sf(coords = c("x", "y"))
  message(paste("\nThere were ", nrow(sampled.points) - nrow(sampled.points.unique),
                "points that were sampled twice. This indicates undersampling of low density regions or oversampling of the border region.\nThis occures as the probability of beeing close to the same points twice is lower in high denisty regions."))

  if (!is.null(n.samples)){
    sample.indexes <- floor(seq(1, nrow(sampled.points.unique),
                                length.out = min(n.samples, nrow(sampled.points.unique))))
    selected.sampled.points <- sampled.points.unique[sample.indexes, ]
    return(selected.sampled.points)
  }

  return(sampled.points.unique)
}

