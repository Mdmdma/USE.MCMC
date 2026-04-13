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

test_that("paSampling rejects NULL grid.res", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSampling(env.rast = r, pres = pres, grid.res = NULL), "'grid.res' must be provided")
})

test_that("paSampling rejects thres out of range", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSampling(env.rast = r, pres = pres, grid.res = 5, thres = 2.0),
               "must be between 0 and 1")
})

test_that("paSampling rejects non-logical sub.ts", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSampling(env.rast = r, pres = pres, grid.res = 5, sub.ts = "yes"),
               "'sub.ts' must be a single logical")
})

test_that("paSampling rejects non-numeric n.tr", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSampling(env.rast = r, pres = pres, grid.res = 5, n.tr = "five"),
               "'n.tr' must be a positive number")
})

# --- Positive tests ---

test_that("paSampling returns sf object", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- paSampling(env.rast = r, pres = pres, grid.res = 3, n.tr = 2,
                       plot_proc = FALSE, verbose = FALSE)
  expect_true(inherits(result, "sf"))
})

test_that("paSampling with sub.ts returns list with obs.tr and obs.ts", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- paSampling(env.rast = r, pres = pres, grid.res = 3, n.tr = 2,
                       sub.ts = TRUE, n.ts = 1, plot_proc = FALSE, verbose = FALSE)
  expect_type(result, "list")
  expect_true("obs.tr" %in% names(result))
  expect_true("obs.ts" %in% names(result))
})
