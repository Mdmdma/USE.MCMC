skip_on_cran()

# --- Input validation tests ---

test_that("thresh.inspect rejects non-raster env.rast", {
  expect_error(thresh.inspect("bad"), "must be a SpatRaster or BasicRaster")
})

test_that("thresh.inspect rejects NULL pres", {
  r <- make_test_raster()
  expect_error(thresh.inspect(r, pres = NULL), "'pres' must be provided")
})

test_that("thresh.inspect rejects non-spatial pres", {
  r <- make_test_raster()
  expect_error(thresh.inspect(r, pres = data.frame(x = 1)), "must be a spatial object")
})

test_that("thresh.inspect rejects thres out of range", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(thresh.inspect(r, pres = pres, thres = 1.5), "must be between 0 and 1")
})

test_that("thresh.inspect rejects negative thres", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 20)
  expect_error(thresh.inspect(r, pres = pres, thres = -0.1), "must be between 0 and 1")
})

# --- Positive tests ---

test_that("thresh.inspect returns ggplot with single threshold", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- thresh.inspect(r, pres = pres, thres = 0.75)
  expect_true(inherits(result, "ggplot"))
})

test_that("thresh.inspect works with vector of thresholds", {
  r <- make_test_raster()
  pres <- make_test_presence_sf(r, 30)
  result <- thresh.inspect(r, pres = pres, thres = c(0.5, 0.75))
  expect_true(inherits(result, "ggplot"))
})
