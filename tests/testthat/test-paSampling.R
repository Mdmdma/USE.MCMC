skip_on_cran()

# --- Input validation tests ---

test_that("paSampling rejects non-raster env.rast", {
  expect_error(paSampling(env.rast = "bad"), "must be a SpatRaster or BasicRaster")
})

test_that("paSampling rejects NULL env.rast", {
  expect_error(paSampling(env.rast = NULL), "'env.rast' must be provided")
})

test_that("paSampling rejects NULL pres", {
  r <- make_test_raster()
  expect_error(paSampling(env.rast = r, pres = NULL), "'pres' must be provided")
})

test_that("paSampling rejects non-spatial pres", {
  r <- make_test_raster()
  expect_error(paSampling(env.rast = r, pres = data.frame(x = 1)), "must be a spatial object")
})

test_that("paSampling rejects thres out of range", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSampling(env.rast = r, pres = pres, grid.res = 5, thres = 2.0),
               "must be between 0 and 1")
})

test_that("paSampling rejects non-numeric n.tr", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSampling(env.rast = r, pres = pres, grid.res = 5, n.tr = "five"),
               "'n.tr' must be a positive number")
})

# --- Cross-sampler argument guidance ---

test_that("paSampling guides parameters that belong to another sampler", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  # a paSamplingNn/Mcmc-only parameter
  expect_error(paSampling(env.rast = r, pres = pres, dimensions = c("PC1", "PC2")),
               "paSamplingNn")
  # a paSamplingMcmc-only parameter
  expect_error(paSampling(env.rast = r, pres = pres, chain.length = 100),
               "paSamplingMcmc")
  # a genuine typo
  expect_error(paSampling(env.rast = r, pres = pres, notarg = 1),
               "unknown argument")
})

# --- Return shape (unified style) ---

test_that("paSampling returns a single sf in the unified style", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- paSampling(env.rast = r, pres = pres, grid.res = 3, n.tr = 2,
                       plot_proc = FALSE, verbose = FALSE)
  expect_true(inherits(result, "sf"))
  # geographic-point geometry, CRS 4326
  expect_equal(sf::st_crs(result)$epsg, 4326L)
  # PC scores are now attribute COLUMNS (not the geometry)
  expect_true(all(c("PC1", "PC2") %in% names(result)))
})

test_that("paSampling works with defaults (no grid.res / shared easy call)", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- paSampling(env.rast = r, pres = pres)
  expect_true(inherits(result, "sf"))
})
