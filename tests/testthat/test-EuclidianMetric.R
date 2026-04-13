# Note: despite the name, euclidianMetric computes Manhattan distance (L1 norm)

test_that("euclidianMetric returns 0 for identical points", {
  a <- data.frame(x = 3, y = 4)
  b <- data.frame(x = 3, y = 4)
  result <- USE.MCMC:::euclidianMetric(a, b, dim = c("x", "y"))
  expect_equal(result, 0)
})

test_that("euclidianMetric computes Manhattan distance correctly", {
  a <- data.frame(x = 0, y = 0)
  b <- data.frame(x = 3, y = 4)
  # Manhattan distance = |3-0| + |4-0| = 7
  result <- USE.MCMC:::euclidianMetric(a, b, dim = c("x", "y"))
  expect_equal(result, 7)
})

test_that("euclidianMetric works with single dimension", {
  a <- data.frame(x = 10, y = 99)
  b <- data.frame(x = 3, y = 99)
  result <- USE.MCMC:::euclidianMetric(a, b, dim = c("x"))
  expect_equal(result, 7)
})

test_that("euclidianMetric works with sf objects", {
  a <- sf::st_as_sf(data.frame(x = 0, y = 0, lon = 1, lat = 1), coords = c("lon", "lat"))
  b <- sf::st_as_sf(data.frame(x = 3, y = 4, lon = 2, lat = 2), coords = c("lon", "lat"))
  result <- USE.MCMC:::euclidianMetric(a, b, dim = c("x", "y"))
  expect_equal(result, 7)
})

test_that("euclidianMetric handles negative values", {
  a <- data.frame(x = -3, y = 2)
  b <- data.frame(x = 3, y = -2)
  # |(-3)-3| + |2-(-2)| = 6 + 4 = 10
  result <- USE.MCMC:::euclidianMetric(a, b, dim = c("x", "y"))
  expect_equal(result, 10)
})
