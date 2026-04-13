# --- Input validation tests ---

test_that("plotDensityLines rejects density=TRUE with NULL dataset", {
  expect_error(
    plotDensityLines(dataset = NULL, density = TRUE),
    "at least one line of data"
  )
})

test_that("plotDensityLines rejects lines=TRUE with too few cols", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    plotDensityLines(dataset = df, lines = TRUE, cols = "x"),
    "At least two columns"
  )
})

test_that("plotDensityLines rejects non-function densityFunction", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    plotDensityLines(dataset = df, density = TRUE, cols = c("x", "y"),
                     densityFunction = "not a function"),
    "density function must be provided"
  )
})

test_that("plotDensityLines rejects bad xlim", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(plotDensityLines(dataset = df, xlim = "bad"), "'xlim' must be a numeric vector")
})

test_that("plotDensityLines rejects bad ylim length", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(plotDensityLines(dataset = df, ylim = c(1, 2, 3)), "'ylim' must be a numeric vector of length 2")
})

test_that("plotDensityLines rejects non-character title", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(plotDensityLines(dataset = df, title = 123), "'title' must be a single character string")
})

test_that("plotDensityLines rejects non-logical lines", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(plotDensityLines(dataset = df, lines = "yes"), "'lines' must be a single logical")
})

test_that("plotDensityLines rejects non-logical density", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(plotDensityLines(dataset = df, density = 1), "'density' must be a single logical")
})

test_that("plotDensityLines rejects non-positive resolution", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(plotDensityLines(dataset = df, resolution = 0), "'resolution' must be a positive number")
})

test_that("plotDensityLines rejects cols not in dataset", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(plotDensityLines(dataset = df, cols = c("x", "z")), "columns not found in 'dataset'")
})

# --- Positive tests ---

test_that("plotDensityLines returns ggplot with default params", {
  df <- data.frame(x = 1:10, y = 1:10)
  result <- plotDensityLines(dataset = df)
  expect_true(inherits(result, "ggplot"))
})

test_that("plotDensityLines returns ggplot with lines=TRUE", {
  df <- data.frame(x = 1:10, y = 1:10)
  result <- plotDensityLines(dataset = df, lines = TRUE, cols = c("x", "y"),
                             xlim = c(0, 12), ylim = c(0, 12))
  expect_true(inherits(result, "ggplot"))
})

test_that("plotDensityLines returns ggplot with density=TRUE", {
  df <- data.frame(x = 1:10, y = 1:10)
  result <- plotDensityLines(dataset = df, density = TRUE, cols = c("x", "y"),
                             densityFunction = alwaysOne, resolution = 3,
                             xlim = c(0, 12), ylim = c(0, 12))
  expect_true(inherits(result, "ggplot"))
})

test_that("plotDensityLines minimal=TRUE returns ggplot", {
  df <- data.frame(x = 1:10, y = 1:10)
  result <- plotDensityLines(dataset = df, minimal = TRUE)
  expect_true(inherits(result, "ggplot"))
})
