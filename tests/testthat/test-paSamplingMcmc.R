skip_on_cran()

# --- Input validation tests ---

test_that("paSamplingMcmc rejects NULL env.data.raster", {
  expect_error(paSamplingMcmc(env.data.raster = NULL), "'env.data.raster' must be provided")
})

test_that("paSamplingMcmc rejects non-raster env.data.raster", {
  expect_error(paSamplingMcmc(env.data.raster = "bad"), "must be a SpatRaster or BasicRaster")
})

test_that("paSamplingMcmc rejects NULL pres", {
  r <- make_test_raster()
  expect_error(paSamplingMcmc(env.data.raster = r, pres = NULL), "'pres' must be provided")
})

test_that("paSamplingMcmc rejects non-spatial pres", {
  r <- make_test_raster()
  expect_error(paSamplingMcmc(env.data.raster = r, pres = data.frame(x = 1)),
               "must be a spatial object")
})

test_that("paSamplingMcmc rejects negative n.samples", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSamplingMcmc(env.data.raster = r, pres = pres, n.samples = -1),
               "'n.samples' must be a positive number")
})

test_that("paSamplingMcmc rejects non-character dimensions", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSamplingMcmc(env.data.raster = r, pres = pres, dimensions = c(1, 2)),
               "'dimensions' must be a character vector")
})

test_that("paSamplingMcmc rejects negative burnIn", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSamplingMcmc(env.data.raster = r, pres = pres, burnIn = -1),
               "'burnIn' must be a non-negative number")
})

test_that("paSamplingMcmc rejects non-logical verbose", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSamplingMcmc(env.data.raster = r, pres = pres, verbose = "yes"),
               "'verbose' must be a single logical")
})

test_that("paSamplingMcmc rejects invalid precomputed.pca", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSamplingMcmc(env.data.raster = r, pres = pres, precomputed.pca = list(a = 1)),
               "must be a list with a '\\$PCs' element")
})

test_that("paSamplingMcmc rejects out-of-range percentile", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(paSamplingMcmc(env.data.raster = r, pres = pres,
                              environmental.cutof.percentile = 2),
               "must be between 0 and 1")
})

# --- Positive tests ---

test_that("paSamplingMcmc returns sf object with small params", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- paSamplingMcmc(
    env.data.raster = r,
    pres = pres,
    n.samples = 10,
    chain.length = 50,
    burnIn = 0,
    num.chains = 1,
    num.cores = 1,
    verbose = FALSE,
    plot_proc = FALSE,
    covariance.correction = 100
  )
  expect_true(inherits(result, "sf"))
  expect_true(nrow(result) <= 10)
})

test_that("paSamplingMcmc with engine = 'cpp' returns sf with density column", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- paSamplingMcmc(
    env.data.raster = r,
    pres = pres,
    n.samples = 50,
    chain.length = 1000,
    burnIn = 200,
    num.chains = 1,
    num.cores = 1,
    verbose = FALSE,
    plot_proc = FALSE,
    covariance.correction = 1,
    engine = "cpp"
  )
  expect_true(inherits(result, "sf"))
  expect_true(nrow(result) <= 50)
  expect_true("density" %in% names(result))
})
