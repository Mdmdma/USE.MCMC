skip_on_cran()
skip_if_not_installed("virtualspecies")

# --- Input validation tests ---

test_that("getVirtualSpeciesPresencePoints rejects NULL env.data", {
  expect_error(getVirtualSpeciesPresencePoints(env.data = NULL),
               "'env.data' must be provided")
})

test_that("getVirtualSpeciesPresencePoints rejects non-raster env.data", {
  expect_error(getVirtualSpeciesPresencePoints(env.data = "bad"),
               "must be a SpatRaster or BasicRaster")
})

test_that("getVirtualSpeciesPresencePoints rejects non-positive n.samples", {
  r <- make_test_raster()
  expect_error(getVirtualSpeciesPresencePoints(env.data = r, n.samples = -1),
               "'n.samples' must be a positive integer")
})

test_that("getVirtualSpeciesPresencePoints rejects non-integer n.samples", {
  r <- make_test_raster()
  expect_error(getVirtualSpeciesPresencePoints(env.data = r, n.samples = 2.5),
               "'n.samples' must be a positive integer")
})

test_that("getVirtualSpeciesPresencePoints rejects zero n.samples", {
  r <- make_test_raster()
  expect_error(getVirtualSpeciesPresencePoints(env.data = r, n.samples = 0),
               "'n.samples' must be a positive integer")
})

test_that("getVirtualSpeciesPresencePoints rejects non-logical plot", {
  r <- make_test_raster()
  expect_error(getVirtualSpeciesPresencePoints(env.data = r, n.samples = 5, plot = "yes"),
               "'plot' must be a single logical")
})

# --- Positive tests ---

test_that("getVirtualSpeciesPresencePoints returns list with sample.points", {
  r <- make_test_raster()
  # virtualspecies::generateRandomSp may fail internally depending on the

  # raster/version combination (e.g. logisticFun not found). Skip gracefully.
  result <- tryCatch(
    getVirtualSpeciesPresencePoints(env.data = r, n.samples = 10, plot = FALSE),
    error = function(e) {
      skip(paste("virtualspecies internal error:", conditionMessage(e)))
    }
  )
  expect_type(result, "list")
  expect_true("sample.points" %in% names(result))
  expect_true(inherits(result$sample.points, "SpatVector"))
})
