skip_on_cran()

# --- Input validation tests ---

test_that("paSamplingNn rejects non-raster env.rast", {
  expect_error(paSamplingNn(env.rast = "bad"), "must be a SpatRaster or BasicRaster")
})

test_that("paSamplingNn rejects NULL env.rast", {
  expect_error(paSamplingNn(env.rast = NULL), "'env.rast' must be provided")
})

test_that("paSamplingNn rejects non-spatial pres", {
  r <- make_test_raster()
  expect_error(paSamplingNn(env.rast = r, pres = data.frame(x = 1)),
               "must be a spatial object")
})

test_that("paSamplingNn rejects non-character dimensions", {
  r <- make_test_raster()
  expect_error(paSamplingNn(env.rast = r, dimensions = c(1, 2)),
               "'dimensions' must be a character vector")
})

test_that("paSamplingNn rejects thres out of range", {
  r <- make_test_raster()
  expect_error(paSamplingNn(env.rast = r, thres = 1.5), "must be between 0 and 1")
})

test_that("paSamplingNn rejects non-logical nn.based.presence.exclusion", {
  r <- make_test_raster()
  expect_error(paSamplingNn(env.rast = r, nn.based.presence.exclusion = "yes"),
               "'nn.based.presence.exclusion' must be a single logical")
})

test_that("paSamplingNn rejects invalid precomputed.pca", {
  r <- make_test_raster()
  expect_error(paSamplingNn(env.rast = r, precomputed.pca = list(a = 1)),
               "must be a list with a '\\$PCs' element")
})

test_that("paSamplingNn rejects non-positive n.samples", {
  r <- make_test_raster()
  expect_error(paSamplingNn(env.rast = r, n.samples = -5),
               "'n.samples' must be NULL or a positive number")
})

# --- Regression: nn-based presence exclusion (previously broken by undefined symbols) ---

test_that("paSamplingNn nn.based.presence.exclusion=TRUE returns an sf object", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- paSamplingNn(env.rast = r, pres = pres, grid.res = 5, n.tr = 2,
                         nn.based.presence.exclusion = TRUE,
                         data.based.distance.threshold = FALSE,
                         n.samples = 10, verbose = FALSE, plot_proc = FALSE)
  expect_true(inherits(result, "sf"))
  expect_true(nrow(result) <= 10)
})

# --- Positive tests ---

test_that("paSamplingNn returns sf object without presence exclusion", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- paSamplingNn(
    env.rast = r, pres = pres, grid.res = 5, n.tr = 2,
    nn.based.presence.exclusion = FALSE,
    data.based.distance.threshold = FALSE,
    n.samples = 20, verbose = FALSE, plot_proc = FALSE
  )
  expect_true(inherits(result, "sf"))
  expect_true(nrow(result) <= 20)
})

test_that("paSamplingNn works without pres (uniform sampling)", {
  r <- make_test_raster()
  result <- paSamplingNn(
    env.rast = r, pres = NULL, grid.res = 5, n.tr = 2,
    nn.based.presence.exclusion = FALSE,
    data.based.distance.threshold = FALSE,
    n.samples = 20, verbose = FALSE, plot_proc = FALSE
  )
  expect_true(inherits(result, "sf"))
})

# --- Arbitrary-dimension support ---

test_that("paSamplingNn rejects fewer than two dimensions", {
  r <- make_test_raster()
  expect_error(paSamplingNn(env.rast = r, dimensions = "PC1"),
               "at least 2 elements")
})

test_that("paSamplingNn samples in 3D with nn-based exclusion + data threshold", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 40)
  result <- paSamplingNn(
    env.rast = r, pres = pres, grid.res = 6, n.tr = 3,
    dimensions = c("PC1", "PC2", "PC3"),
    nn.based.presence.exclusion = TRUE,
    data.based.distance.threshold = TRUE,
    n.samples = 10, verbose = FALSE, plot_proc = FALSE
  )
  expect_true(inherits(result, "sf"))
  expect_true(nrow(result) <= 10)
  # geometry stays geographic (2D) regardless of the PC-space dimensionality
  expect_true(all(c("PC1", "PC2", "PC3") %in% names(result)))
})

test_that("paSamplingNn samples in 4D (compute-bound, more candidates)", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 40)
  result <- paSamplingNn(
    env.rast = r, pres = pres,
    dimensions = c("PC1", "PC2", "PC3", "PC4"),
    nn.based.presence.exclusion = TRUE,
    data.based.distance.threshold = TRUE,
    n.candidates = 4000, n.samples = 8, verbose = FALSE, plot_proc = FALSE
  )
  expect_true(inherits(result, "sf"))
  expect_true(nrow(result) <= 8)
})
