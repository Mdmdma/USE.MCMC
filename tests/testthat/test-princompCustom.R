# --- Input validation tests ---

test_that("princompCustom rejects invalid covmat list", {
  expect_error(
    princompCustom(covmat = list(wrong = 1)),
    "not a valid covariance list"
  )
})

test_that("princompCustom rejects non-numeric data", {
  expect_error(
    princompCustom(covmat = matrix(c("a", "b", "c", "d"), 2, 2)),
    "numerical variables"
  )
})

test_that("princompCustom rejects cor=TRUE with constant variable", {
  # Create cov matrix with zero variance for one variable
  cv <- matrix(c(1, 0, 0, 0), 2, 2)
  expect_error(
    princompCustom(covmat = cv, cor = TRUE),
    "constant variable"
  )
})

test_that("princompCustom warns when both x and covmat are supplied", {
  x <- matrix(rnorm(50), 25, 2)
  cv <- cov(x)
  expect_warning(
    princompCustom(x = x, covmat = cv),
    "both.*x.*and.*covmat"
  )
})

# --- Positive tests ---

test_that("princompCustom returns princomp object from covariance matrix", {
  set.seed(42)
  x <- matrix(rnorm(200), 100, 2)
  cv <- cov(x)
  result <- princompCustom(covmat = cv)
  expect_s3_class(result, "princomp")
  expect_true(!is.null(result$sdev))
  expect_true(!is.null(result$loadings))
})

test_that("princompCustom returns correct number of components", {
  set.seed(42)
  x <- matrix(rnorm(300), 100, 3)
  result <- princompCustom(x)
  expect_length(result$sdev, 3)
  expect_equal(ncol(result$loadings), 3)
})

test_that("princompCustom sdev values are non-negative", {
  set.seed(42)
  x <- matrix(rnorm(200), 100, 2)
  result <- princompCustom(x)
  expect_true(all(result$sdev >= 0))
})

test_that("fix_sign=TRUE makes first element of loadings non-negative", {
  set.seed(42)
  x <- matrix(rnorm(300), 100, 3)
  result <- princompCustom(x, fix_sign = TRUE)
  first_elements <- result$loadings[1, ]
  expect_true(all(first_elements >= 0))
})

test_that("princompCustom computes scores when scores=TRUE", {
  set.seed(42)
  x <- matrix(rnorm(200), 100, 2)
  result <- princompCustom(x, scores = TRUE)
  expect_equal(nrow(result$scores), 100)
  expect_equal(ncol(result$scores), 2)
})

test_that("princompCustom works with cov.wt list", {
  set.seed(42)
  x <- matrix(rnorm(200), 100, 2)
  covmat <- cov.wt(x)
  result <- princompCustom(covmat = covmat)
  expect_s3_class(result, "princomp")
})
