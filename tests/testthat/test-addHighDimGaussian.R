# --- Input validation tests ---

test_that("addHighDimGaussian rejects dim = 0", {
  expect_error(addHighDimGaussian(dim = 0), "'dim' must be a positive integer")
})

test_that("addHighDimGaussian rejects negative dim", {
  expect_error(addHighDimGaussian(dim = -1), "'dim' must be a positive integer")
})

test_that("addHighDimGaussian rejects non-integer dim", {
  expect_error(addHighDimGaussian(dim = 1.5), "'dim' must be a positive integer")
})

test_that("addHighDimGaussian rejects non-numeric dim", {
  expect_error(addHighDimGaussian(dim = "a"), "'dim' must be a positive integer")
})

test_that("addHighDimGaussian rejects mismatched mean.vec length", {
  expect_error(
    addHighDimGaussian(dim = 2, mean.vec = c(0, 0, 0)),
    "'mean.vec' length must equal 'dim'"
  )
})

test_that("addHighDimGaussian rejects non-numeric mean.vec", {
  expect_error(
    addHighDimGaussian(dim = 2, mean.vec = c("a", "b")),
    "'mean.vec' must be numeric"
  )
})

test_that("addHighDimGaussian rejects mismatched cov.mat dimensions", {
  expect_error(
    addHighDimGaussian(dim = 2, cov.mat = diag(3)),
    "'cov.mat' must be a 2x2 matrix, got 3x3"
  )
})

test_that("addHighDimGaussian rejects non-matrix cov.mat", {
  expect_error(
    addHighDimGaussian(dim = 2, cov.mat = c(1, 0, 0, 1)),
    "'cov.mat' must be a numeric matrix"
  )
})

# --- Positive tests ---

test_that("addHighDimGaussian returns a function", {
  fn <- addHighDimGaussian(dim = 2)
  expect_true(is.function(fn))
})

test_that("returned function adds noise to a numeric vector", {
  set.seed(42)
  fn <- addHighDimGaussian(dim = 2)
  point <- c(0, 0)
  result <- fn(point, covariance.adjuster = 1, dim = c("PC1", "PC2"))
  # Should have modified the values (almost certainly not exactly 0)
  expect_true(result[1] != 0 || result[2] != 0)
  expect_length(result, 2)
})

test_that("returned function works in 1D", {
  set.seed(42)
  fn <- addHighDimGaussian(dim = 1)
  point <- 5
  result <- fn(point, covariance.adjuster = 1, dim = c("x"))
  expect_true(is.numeric(result))
  expect_length(result, 1)
})

test_that("covariance.adjuster scales the noise", {
  fn <- addHighDimGaussian(dim = 1, cov.mat = matrix(1))
  # Large adjuster -> larger noise variance
  set.seed(42)
  results_small <- replicate(100, {
    fn(0, covariance.adjuster = 0.001, dim = "x")
  })
  set.seed(42)
  results_large <- replicate(100, {
    fn(0, covariance.adjuster = 100, dim = "x")
  })
  expect_true(sd(results_large) > sd(results_small))
})
