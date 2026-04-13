# gaussianMixtureDensityFunction is an internal function

test_that("gaussianMixtureDensityFunction returns a function", {
  set.seed(42)
  df <- data.frame(PC1 = rnorm(100), PC2 = rnorm(100))
  sf_data <- sf::st_as_sf(df, coords = c("PC1", "PC2"))
  # Need columns accessible after st_drop_geometry, so add them as regular cols too
  sf_data$PC1 <- df$PC1
  sf_data$PC2 <- df$PC2
  fn <- USE.MCMC:::gaussianMixtureDensityFunction(sf_data, dim = c("PC1", "PC2"))
  expect_true(is.function(fn))
})

test_that("returned function returns numeric scalar", {
  set.seed(42)
  df <- data.frame(PC1 = rnorm(100), PC2 = rnorm(100))
  sf_data <- sf::st_as_sf(df, coords = c("PC1", "PC2"))
  sf_data$PC1 <- df$PC1
  sf_data$PC2 <- df$PC2
  fn <- USE.MCMC:::gaussianMixtureDensityFunction(sf_data, dim = c("PC1", "PC2"))
  point <- c(0, 0)
  result <- fn(point)
  expect_true(is.numeric(result))
  expect_length(result, 1)
})

test_that("returned function returns threshold/1000 for far-out points", {
  set.seed(42)
  df <- data.frame(PC1 = rnorm(100), PC2 = rnorm(100))
  sf_data <- sf::st_as_sf(df, coords = c("PC1", "PC2"))
  sf_data$PC1 <- df$PC1
  sf_data$PC2 <- df$PC2
  threshold <- 0.01
  fn <- USE.MCMC:::gaussianMixtureDensityFunction(sf_data, dim = c("PC1", "PC2"), threshold = threshold)
  far_point <- c(1000, 1000)
  result <- fn(far_point)
  expect_equal(result, threshold / 1000)
})

test_that("returned function returns 1 for points within distribution", {
  set.seed(42)
  df <- data.frame(PC1 = rnorm(200), PC2 = rnorm(200))
  sf_data <- sf::st_as_sf(df, coords = c("PC1", "PC2"))
  sf_data$PC1 <- df$PC1
  sf_data$PC2 <- df$PC2
  fn <- USE.MCMC:::gaussianMixtureDensityFunction(sf_data, dim = c("PC1", "PC2"), threshold = 0.001)
  center_point <- c(0, 0)
  result <- fn(center_point)
  expect_equal(result, 1)
})
