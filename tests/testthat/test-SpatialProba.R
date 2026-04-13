# --- Input validation tests ---

test_that("SpatialProba rejects NULL coefs", {
  r <- make_test_raster()
  expect_error(SpatialProba(coefs = NULL, env.rast = r), "'coefs' must be provided")
})

test_that("SpatialProba rejects non-numeric coefs", {
  r <- make_test_raster()
  expect_error(SpatialProba(coefs = c(a = "x"), env.rast = r), "'coefs' must be a numeric vector")
})

test_that("SpatialProba rejects unnamed coefs", {
  r <- make_test_raster()
  expect_error(SpatialProba(coefs = c(1, 2), env.rast = r), "named vector")
})

test_that("SpatialProba rejects NULL env.rast", {
  coefs <- c(intercept = 0.1, bio1 = 0.5)
  expect_error(SpatialProba(coefs = coefs, env.rast = NULL), "'env.rast' must be provided")
})

test_that("SpatialProba rejects non-raster env.rast", {
  coefs <- c(intercept = 0.1, bio1 = 0.5)
  expect_error(SpatialProba(coefs = coefs, env.rast = "bad"), "must be a SpatRaster or BasicRaster")
})

test_that("SpatialProba rejects non-logical marginalPlots", {
  r <- make_test_raster()
  layer_name <- names(r)[1]
  coefs <- c(intercept = 0.1)
  coefs[layer_name] <- 0.5
  expect_error(SpatialProba(coefs = coefs, env.rast = r, marginalPlots = "yes"),
               "'marginalPlots' must be a single logical")
})

test_that("SpatialProba rejects coefs names not in raster", {
  r <- make_test_raster()
  coefs <- c(intercept = 0.1, nonexistent_var = 0.5)
  expect_error(SpatialProba(coefs = coefs, env.rast = r), "not all names")
})

test_that("SpatialProba rejects quadr_term without quad prefix in coefs", {
  r <- make_test_raster()
  layer_name <- names(r)[1]
  coefs <- c(intercept = 0.1)
  coefs[layer_name] <- 0.5
  expect_error(
    SpatialProba(coefs = coefs, env.rast = r, quadr_term = layer_name),
    "no term with prefix"
  )
})

# --- Positive tests ---

test_that("SpatialProba returns SpatRaster with marginalPlots=FALSE", {
  r <- make_test_raster()
  layer_name <- names(r)[1]
  coefs <- c(intercept = 0.1)
  coefs[layer_name] <- 0.01
  result <- SpatialProba(coefs = coefs, env.rast = r, marginalPlots = FALSE)
  expect_true(inherits(result, "SpatRaster"))
  expect_equal(names(result), "TrueProba")
})

test_that("SpatialProba probability values are between 0 and 1", {
  r <- make_test_raster()
  layer_name <- names(r)[1]
  coefs <- c(intercept = 0.1)
  coefs[layer_name] <- 0.01
  result <- SpatialProba(coefs = coefs, env.rast = r, marginalPlots = FALSE)
  vals <- terra::values(result, na.rm = TRUE)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("SpatialProba marginalPlots=TRUE errors when quadr_term is NULL (known bug)", {
  # Pre-existing bug: marginal effects code path at line 88-106 unconditionally

  # accesses quadr_term columns, failing with "subscript out of bounds" when NULL
  r <- make_test_raster()
  layer_names <- names(r)[1:2]
  coefs <- c(intercept = 0.1)
  coefs[layer_names[1]] <- 0.01
  coefs[layer_names[2]] <- -0.005
  expect_error(SpatialProba(coefs = coefs, env.rast = r, marginalPlots = TRUE),
               "subscript out of bounds")
})

test_that("SpatialProba returns list with marginalPlots=TRUE and quadr_term", {
  r <- make_test_raster()
  layer_names <- names(r)[1:2]
  coefs <- c(intercept = 0.1)
  coefs[layer_names[1]] <- 0.01
  coefs[paste0("quad_", layer_names[1])] <- -0.001
  coefs[layer_names[2]] <- -0.005
  result <- SpatialProba(coefs = coefs, env.rast = r,
                         quadr_term = layer_names[1], marginalPlots = TRUE)
  expect_type(result, "list")
  expect_true("rast" %in% names(result))
  expect_true("margEff" %in% names(result))
  expect_true(inherits(result$rast, "SpatRaster"))
})
