# --- Input validation tests ---

test_that("optimalDistanceThresholdNn rejects NULL env.data", {
  expect_error(optimalDistanceThresholdNn(env.data = NULL), "'env.data' must be provided")
})

test_that("optimalDistanceThresholdNn rejects non-dataframe env.data", {
  expect_error(optimalDistanceThresholdNn(env.data = "bad"), "'env.data' must be a data.frame")
})

test_that("optimalDistanceThresholdNn rejects dimensions not in data", {
  df <- data.frame(PC1 = 1:10, PC2 = 1:10)
  expect_error(
    optimalDistanceThresholdNn(env.data = df, dimensions = c("PC1", "PC3")),
    "columns not found in 'env.data'"
  )
})

test_that("optimalDistanceThresholdNn rejects non-integer index.for.cutof", {
  df <- data.frame(PC1 = 1:10, PC2 = 1:10)
  expect_error(
    optimalDistanceThresholdNn(env.data = df, index.for.cutof = 1.5),
    "'index.for.cutof' must be a positive integer"
  )
})

test_that("optimalDistanceThresholdNn rejects non-integer num.neighbors", {
  df <- data.frame(PC1 = 1:10, PC2 = 1:10)
  expect_error(
    optimalDistanceThresholdNn(env.data = df, num.neighbors = 0),
    "'num.neighbors' must be a positive integer"
  )
})

# --- Positive tests ---

test_that("optimalDistanceThresholdNn returns single positive numeric", {
  set.seed(42)
  df <- data.frame(PC1 = rnorm(50), PC2 = rnorm(50))
  result <- optimalDistanceThresholdNn(env.data = df, dimensions = c("PC1", "PC2"))
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(result > 0)
})

test_that("optimalDistanceThresholdNn works with sf input", {
  pc_sf <- make_test_pc_sf(50)
  pc_sf$PC1 <- sf::st_coordinates(pc_sf)[, 1]
  pc_sf$PC2 <- sf::st_coordinates(pc_sf)[, 2]
  # Use x and y as stand-in dimensions (they exist as regular columns in the sf)
  set.seed(42)
  df <- sf::st_drop_geometry(pc_sf)
  result <- optimalDistanceThresholdNn(env.data = df, dimensions = c("PC1", "PC2"))
  expect_true(is.numeric(result))
  expect_true(result > 0)
})
