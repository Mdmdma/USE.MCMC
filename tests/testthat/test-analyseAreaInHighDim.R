skip_on_cran()

# --- Input validation tests ---

test_that("analyseAreaInHigherDim rejects NULL data.df", {
  expect_error(USE.MCMC:::analyseAreaInHigherDim(data.df = NULL),
               "'data.df' must be a data.frame")
})

test_that("analyseAreaInHigherDim rejects non-data.frame", {
  expect_error(USE.MCMC:::analyseAreaInHigherDim(data.df = "bad"),
               "'data.df' must be a data.frame")
})

test_that("analyseAreaInHigherDim rejects empty data.frame", {
  expect_error(USE.MCMC:::analyseAreaInHigherDim(data.df = data.frame()),
               "must have at least one row")
})

test_that("analyseAreaInHigherDim rejects non-character dim.to.analyse", {
  df <- data.frame(PC1 = 1:10, PC2 = 1:10, PC3 = 1:10)
  expect_error(USE.MCMC:::analyseAreaInHigherDim(df, dim.to.analyse = 1),
               "'dim.to.analyse' must be a character vector")
})

test_that("analyseAreaInHigherDim rejects non-character dim.to.compress", {
  df <- data.frame(PC1 = 1:10, PC2 = 1:10, PC3 = 1:10)
  expect_error(USE.MCMC:::analyseAreaInHigherDim(df, dim.to.compress = c(1, 2)),
               "'dim.to.compress' must be a character vector")
})

test_that("analyseAreaInHigherDim rejects missing columns", {
  df <- data.frame(PC1 = 1:10, PC2 = 1:10)
  expect_error(USE.MCMC:::analyseAreaInHigherDim(df, dim.to.analyse = c("PC3")),
               "Columns not found in 'data.df': PC3")
})

test_that("analyseAreaInHigherDim rejects non-positive num.intervals", {
  df <- data.frame(PC1 = 1:10, PC2 = 1:10, PC3 = 1:10)
  expect_error(USE.MCMC:::analyseAreaInHigherDim(df, num.intervals = 0),
               "'num.intervals' must be a positive integer")
})

test_that("analyseAreaInHigherDim rejects non-integer num.intervals", {
  df <- data.frame(PC1 = 1:10, PC2 = 1:10, PC3 = 1:10)
  expect_error(USE.MCMC:::analyseAreaInHigherDim(df, num.intervals = 2.5),
               "'num.intervals' must be a positive integer")
})

# --- Positive tests ---

test_that("analyseAreaInHigherDim returns data.frame with correct structure", {
  set.seed(42)
  df <- data.frame(
    PC1 = rnorm(200),
    PC2 = rnorm(200),
    PC3 = rnorm(200)
  )
  result <- USE.MCMC:::analyseAreaInHigherDim(df, num.intervals = 5)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5)
  expect_true("interval.average" %in% names(result))
  expect_true("area" %in% names(result))
  # area should be normalized (sums to 1, ignoring intervals with 0 area)
  expect_equal(sum(result$area), 1, tolerance = 1e-10)
})

test_that("analyseAreaInHigherDim works with custom dimensions", {
  set.seed(42)
  df <- data.frame(
    dim_a = rnorm(200),
    dim_b = rnorm(200),
    dim_c = rnorm(200)
  )
  result <- USE.MCMC:::analyseAreaInHigherDim(
    df,
    dim.to.analyse = c("dim_c"),
    dim.to.compress = c("dim_a", "dim_b"),
    num.intervals = 3
  )
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)
})
