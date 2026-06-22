# Standalone validation of the dimension-general NN sampling CORE (FNN + base R only;
# no terra/sf needed). Mirrors the algorithmic block of paSamplingNn() verbatim, run
# on a synthetic d-dimensional cloud, to validate behaviour for d > 2.
suppressMessages(library(FNN))
set.seed(1)
fail <- 0L
ok <- function(cond, msg) { cat(sprintf("[%s] %s\n", ifelse(cond, "PASS", "FAIL"), msg)); if (!cond) fail <<- fail + 1L }

# the real correction helper, sourced from the package file
src <- parse("R/optimalDistanceThresholdNn.R"); eval(src[[length(src)]])

# --- core, copied verbatim from paSamplingNn() body (operating on a matrix) ----------
nn_core <- function(pc.mat, dimensions, grid.res = 10, n.tr = 5, n.samples = NULL,
                    n.candidates = NULL, data.based.distance.threshold = TRUE,
                    dim.correction = "voronoi", thres = 0.75, pres.mat = NULL) {
  d <- length(dimensions)
  axis.min  <- apply(pc.mat, 2, min)
  axis.span <- apply(pc.mat, 2, function(x) diff(range(x)))
  step <- axis.span / max(1, grid.res)
  if (is.null(n.candidates)) n.candidates <- max(as.integer(round(grid.res^2 * n.tr)), 1L)
  n.candidates <- as.integer(round(n.candidates))

  if (data.based.distance.threshold) {
    nnd <- as.vector(FNN::knn.dist(pc.mat, k = 3))
    idx <- 5
    top.cutoff <- -sort(-nnd, partial = idx)[idx]
    distance.threshold <- (top.cutoff / 2) * .nnThresholdDimFactor(d, dim.correction)
  } else {
    distance.threshold <- (max(step) / 2) * .nnThresholdDimFactor(d, dim.correction)
  }

  target <- if (!is.null(n.samples)) n.samples else NA_real_
  max.batches <- if (!is.null(n.samples)) 25L else 1L
  kept.index <- integer(0); kept.dist <- numeric(0); n.accepted <- 0L
  for (b in seq_len(max.batches)) {
    u <- matrix(stats::runif(n.candidates * d), ncol = d)
    u <- sweep(sweep(u, 2, axis.span, `*`), 2, axis.min, `+`)
    nn <- FNN::get.knnx(pc.mat, u, k = 1)
    keep <- nn$nn.dist[, 1] < distance.threshold
    if (any(keep)) {
      n.accepted <- n.accepted + sum(keep)
      kept.index <- c(kept.index, nn$nn.index[keep, 1])
      kept.dist  <- c(kept.dist,  nn$nn.dist[keep, 1])
      first <- !duplicated(kept.index); kept.index <- kept.index[first]; kept.dist <- kept.dist[first]
    }
    if (!is.na(target) && length(kept.index) >= target) break
  }
  if (length(kept.index) == 0L) return(list(index = integer(0), threshold = distance.threshold, n.dup = 0L))

  outside <- rep(TRUE, length(kept.index))
  if (!is.null(pres.mat)) {
    cand.mat <- pc.mat[kept.index, , drop = FALSE]
    k.neighbors <- min(100L, nrow(pres.mat))
    pres.self  <- FNN::get.knnx(pres.mat, pres.mat, k = min(k.neighbors + 1L, nrow(pres.mat)))
    pres.proxy <- rowMeans(pres.self$nn.dist[, -1, drop = FALSE])
    cand.knn   <- FNN::get.knnx(pres.mat, cand.mat, k = k.neighbors)
    cand.proxy <- rowMeans(cand.knn$nn.dist)
    proxy.cut  <- stats::quantile(pres.proxy, 1 - thres)
    outside    <- cand.proxy > proxy.cut
  }
  idx.final <- kept.index[outside]
  if (!is.null(n.samples) && length(idx.final) > 0) {
    si <- floor(seq(1, length(idx.final), length.out = min(n.samples, length(idx.final))))
    idx.final <- idx.final[si]
  }
  list(index = idx.final, threshold = distance.threshold, n.dup = n.accepted - length(kept.index))
}

# ----- a synthetic cloud with strongly NON-uniform density (dense core + sparse halo) -----
make_cloud <- function(d, n = 6000) {
  n1 <- round(n * 0.8); n2 <- n - n1
  core <- matrix(rnorm(n1 * d, 0, 0.5), ncol = d)      # dense core
  halo <- matrix(rnorm(n2 * d, 0, 3.0), ncol = d)      # sparse halo
  rbind(core, halo)
}
local_density <- function(mat, pts) {              # 1 / mean distance to 10 NN in the cloud
  dd <- FNN::get.knnx(mat, pts, k = 10)$nn.dist
  1 / (rowMeans(dd) + 1e-9)
}

for (d in c(2, 3, 5)) {
  dims <- paste0("PC", seq_len(d))
  cl <- make_cloud(d)
  res <- nn_core(cl, dims, grid.res = 12, n.tr = 4, n.samples = 400)
  ok(length(res$index) > 0, sprintf("d=%d: produced %d pseudo-absences", d, length(res$index)))
  ok(all(res$index %in% seq_len(nrow(cl))), sprintf("d=%d: all indices valid", d))
  ok(!any(duplicated(res$index)), sprintf("d=%d: indices unique (deduped)", d))
  ok(length(res$index) <= 400, sprintf("d=%d: respects n.samples cap", d))

  # CORE CLAIM: remap+threshold flattens density vs naive random draws from the cloud.
  set.seed(7)
  rnd <- sample(nrow(cl), length(res$index))
  dens_sel <- local_density(cl, cl[res$index, , drop = FALSE])
  dens_rnd <- local_density(cl, cl[rnd, , drop = FALSE])
  ok(median(dens_sel) < median(dens_rnd),
     sprintf("d=%d: selected points sit in LOWER-density regions than random (median %.3f < %.3f) -> flattening",
             d, median(dens_sel), median(dens_rnd)))

  # threshold correction grows with d
  if (d == 2) thr2 <- res$threshold
}

# threshold correction: same cloud, voronoi vs none at d=5
cl5 <- make_cloud(5); dims5 <- paste0("PC", 1:5)
t_vor <- nn_core(cl5, dims5, n.samples = 50)$threshold
t_non <- nn_core(cl5, dims5, n.samples = 50, dim.correction = "none")$threshold
ok(abs(t_vor / t_non - sqrt(5/2)) < 1e-6, sprintf("d=5: voronoi threshold = none * sqrt(5/2) (%.4f vs %.4f)", t_vor, t_non))

# presence exclusion DIRECTION: presences = a dense blob; candidates in the blob must be dropped
set.seed(3)
cl <- make_cloud(3); dims <- paste0("PC", 1:3)
pres <- matrix(rnorm(300 * 3, 2, 0.3), ncol = 3)        # tight presence blob centred at (2,2,2)
res_noex <- nn_core(cl, dims, n.samples = 1000)
res_ex   <- nn_core(cl, dims, n.samples = 1000, pres.mat = pres)
# fraction of kept points near the presence centre should drop sharply after exclusion
near_centre <- function(idx) mean(sqrt(rowSums((cl[idx, , drop=FALSE] - 2)^2)) < 1)
ok(near_centre(res_ex$index) < near_centre(res_noex$index),
   sprintf("3D exclusion removes presence-core points (near-centre frac %.3f -> %.3f)",
           near_centre(res_noex$index), near_centre(res_ex$index)))

# empty-survivor path: impossibly tight threshold via tiny dim.correction multiplier yields no/ few survivors gracefully
res_tight <- nn_core(make_cloud(4), paste0("PC",1:4), n.samples = 10, dim.correction = 1e-9)
ok(is.integer(res_tight$index), "tiny threshold returns gracefully (no error)")

cat(sprintf("\n==== %s (%d failures) ====\n", ifelse(fail == 0, "ALL CORE CHECKS PASSED", "FAILURES"), fail))
quit(status = ifelse(fail == 0, 0, 1))
