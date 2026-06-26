#' Sampling pseudo-absences for the training and testing datasets.  
#'
#' \code{paSampling} performs a two-step procedure for uniformly sampling pseudo-absences within the environmental space.
#' In the initial step, a kernel-based filter is utilized to determine the subset of the environmental space that will be subsequently sampled. The kernel-based filter calculates the probability function based on the presence observations, enabling the identification of areas within the environmental space that likely exhibit suitable conditions for the species. To achieve this, a probability threshold value is utilized to assign observations to the corresponding portion of the environmental space. These areas, deemed to have suitable environmental conditions, are excluded from the subsequent uniform sampling process conducted in the second step using the \code{uniformSampling} function, which is internally called. 
#' The bandwidth of the kernel can be automatically estimated from the presence observations or directly set by the user, providing flexibility in determining the scope and precision of the filter.
#' 
#' Being designed with species distribution models in mind, \code{paSampling} allows collectively sampling pseudo-absences for both the training and testing dataset (optional). In both cases, the user must provide a number of observations that will be sampled in each cell of the sampling grid (\code{n.tr}: points for the training dataset; \code{n.ts}: points for the testing dataset). 
#' Note that the optimal resolution of the sampling grid can be found using the \code{optimRes} function. Also, note that the number of pseudo-absences eventually sampled in each cell by the internally-called \code{uniformSampling} function depends on the spatial configuration of the observations within the environmental space. Indeed, in most cases some cells of the sampling grid will be empty (i.e., those at the boundary of the environmental space). For this reason, the number of pseudo-absences returned by \code{paSampling} is likely to be lower than the product between the number of cells of the sampling grid and \code{n.tr}(or \code{n.ts}).
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
#' @importFrom stats na.omit quantile
#' @return An sf object with the coordinates of the pseudo-absences both in the geographical and environmental space.
#' @examples
#' \donttest{
#' env <- terra::rast(USE.MCMC::Worldclim_tmp, type = "xyz")
#' df  <- terra::as.data.frame(env, xy = TRUE, na.rm = TRUE)
#' set.seed(1)
#' pres <- sf::st_as_sf(df[sample(nrow(df), 50), ], coords = c("x", "y"), crs = 4326)
#' pa <- paSampling(env.rast = env, pres = pres, grid.res = 10, n.tr = 2)
#' }
#' @export
#'
paSampling <- function (env.rast=NULL, pres = NULL, thres = 0.75, H = NULL, grid.res = NULL,
                         n.tr = 5, sub.ts = FALSE, n.ts = 5, prev = NULL, plot_proc = FALSE,
                         verbose = FALSE)
{
  # Input validation
  check_raster_input(env.rast, "env.rast")
  check_spatial_points(pres, "pres")
  if (is.null(grid.res)) {
    stop("'grid.res' must be provided as a positive integer (got NULL)", call. = FALSE)
  }
  if (!is.numeric(grid.res) || length(grid.res) != 1 || grid.res < 1) {
    stop(paste0("'grid.res' must be a positive number, got ", deparse(grid.res)), call. = FALSE)
  }
  check_in_range(thres, "thres", min_val = 0, max_val = 1)
  if (!is.numeric(n.tr) || length(n.tr) != 1 || n.tr < 1) {
    stop(paste0("'n.tr' must be a positive number, got ", deparse(n.tr)), call. = FALSE)
  }
  if (!is.logical(sub.ts) || length(sub.ts) != 1) {
    stop(paste0("'sub.ts' must be a single logical value, got '", paste(class(sub.ts), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop(paste0("'verbose' must be a single logical value, got '", paste(class(verbose), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.logical(plot_proc) || length(plot_proc) != 1) {
    stop(paste0("'plot_proc' must be a single logical value, got '", paste(class(plot_proc), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.null(prev)) {
    check_in_range(prev, "prev", min_val = 0, max_val = 1)
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
  message("Computing PCA and presences kernel density estimation in the multivariate space")
  rpc <- rastPCA(env.rast, stand = TRUE)
  id_rast <- terra::rast(vals= 1:terra::ncell(env.rast),
                         names ="myID",
                         extent = terra::ext(env.rast),
                         nrows = terra::nrow(env.rast),
                         ncols = terra::ncol(env.rast), 
                         crs = terra::crs(env.rast)
  )
  abio.st <- terra::as.data.frame(c(id_rast, env.rast))
  dt <- terra::as.data.frame(c(id_rast, rpc$PCs[[c("PC1", "PC2")]]), xy = TRUE)
  PC12occ <- terra::extract(id_rast,  occ.vec, cells = FALSE, df = FALSE, ID=FALSE)[,1] 
  PC12ex <- na.omit(data.frame(dt, PA= ifelse(dt$myID %in% PC12occ, 1, 0)))
  if (is.null(H)) {
    H <- ks::Hpi(x = PC12ex[, c("PC1", "PC2")])
  }
  estimate <- data.frame(KDE = ks::kde(PC12ex[PC12ex$PA == 1, c("PC1", "PC2")], 
                                       eval.points = PC12ex[PC12ex$PA ==1, c("PC1", "PC2")], h = H)$estimate, 
                         PC12ex[PC12ex$PA == 1, c("PC1", "PC2", "myID", "PA")])
  quantP <- quantile(estimate[, "KDE"], thres)
  estimate$percP <- ifelse(estimate$KDE <= unname(quantP[1]), "out", "in")
  estimate <- merge(x = PC12ex, y = estimate[estimate$PA == 1,c("myID","percP")], 
                    by = "myID", all.x = TRUE) 
  estimate$percP <- ifelse(is.na(estimate$percP), "pabs",estimate$percP)
  chull <- sf::st_as_sf(subset(estimate, estimate$percP=="in", select=c( "PC1","PC2" )), coords=c( "PC1","PC2" )) 
  chull <- sf::st_union(chull)
  chull <- sf::st_convex_hull(chull)
  estimate <- sf::st_as_sf(estimate, coords=c( "PC1","PC2" ))
  within_chull <- estimate[sf::st_within(estimate, chull, sparse = FALSE), ]
  estimate$percP <- ifelse(estimate$myID %in% within_chull$myID, "in", estimate$percP)
  estimate <- subset(estimate, estimate$percP =="pabs")
  fullDB.sp <- merge(x = estimate, y = abio.st, 
                     by = "myID", all.x = TRUE)
  fullDB.sp <- sf::st_as_sf(fullDB.sp, coords = c("PC1", "PC2"))
  if (is.null(prev)) {
    myPas <- NULL
  }else {
    myPas <- floor(nrow(pres)/prev)
  }
  message("\nPerforming pseudo-absences sampling in the environmental space\n")
  Res <- uniformSampling(sdf = fullDB.sp, grid.res = grid.res, 
                         n.tr = n.tr, n.prev = myPas, sub.ts = sub.ts, n.ts = n.ts, 
                         plot_proc = plot_proc, verbose = verbose)
  if (sub.ts) {
    message("\n", paste(nrow(Res$obs.tr), "training pseudo-absences sampled in the environmental space, \n and", 
                        nrow(Res$obs.ts), "testing pseudo-absences sampled in the environmental space.", 
                        sep = " "), "\n")
    message("\n", paste("Estimated final prevalence of", 
                        round(nrow(pres)/nrow(Res$obs.tr), 2), "instead of", 
                        prev), "\n")
  }
  else {
    message("\n", paste(nrow(Res), "training pseudo-absences sampled in the environmental space", 
                        sep = " "), "\n")
    message("\n", paste("Estimated final prevalence of", 
                        round(nrow(pres)/nrow(Res), 2), "instead of", prev), 
            "\n")
  }
  return(Res)
}
