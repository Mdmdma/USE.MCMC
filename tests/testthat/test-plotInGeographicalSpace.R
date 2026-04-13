skip_if_not_installed("tidyterra")

# --- Input validation tests ---

test_that("plotInGeographicalSpace rejects NULL raster", {
  expect_error(plotInGeographicalSpace(presence.distribution.raster = NULL),
               "'presence.distribution.raster' must be provided")
})

test_that("plotInGeographicalSpace rejects non-raster input", {
  expect_error(plotInGeographicalSpace(presence.distribution.raster = "bad"),
               "must be a SpatRaster or BasicRaster")
})

test_that("plotInGeographicalSpace rejects NULL presence.points", {
  r <- make_test_raster()
  expect_error(plotInGeographicalSpace(presence.distribution.raster = r[[1]],
                                       presence.points = NULL),
               "'presence.points' must be provided")
})

test_that("plotInGeographicalSpace rejects non-sf presence.points", {
  r <- make_test_raster()
  expect_error(plotInGeographicalSpace(presence.distribution.raster = r[[1]],
                                       presence.points = data.frame(x = 1)),
               "'presence.points' must be an sf object")
})

test_that("plotInGeographicalSpace rejects NULL absence.points", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 5)
  expect_error(plotInGeographicalSpace(presence.distribution.raster = r[[1]],
                                       presence.points = pres,
                                       absence.points = NULL),
               "'absence.points' must be provided")
})

test_that("plotInGeographicalSpace rejects non-logical minimal", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 5)
  expect_error(plotInGeographicalSpace(presence.distribution.raster = r[[1]],
                                       presence.points = pres,
                                       absence.points = pres,
                                       minimal = "yes"),
               "'minimal' must be a single logical")
})

# --- Positive tests ---

test_that("plotInGeographicalSpace returns ggplot", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 10)
  abs <- make_test_presence_sf(r, 10)
  result <- plotInGeographicalSpace(
    presence.distribution.raster = r[[1]],
    presence.points = pres,
    absence.points = abs
  )
  expect_true(inherits(result, "ggplot"))
})

test_that("plotInGeographicalSpace works with minimal=TRUE", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 10)
  abs <- make_test_presence_sf(r, 10)
  result <- plotInGeographicalSpace(
    presence.distribution.raster = r[[1]],
    presence.points = pres,
    absence.points = abs,
    minimal = TRUE
  )
  expect_true(inherits(result, "ggplot"))
})
