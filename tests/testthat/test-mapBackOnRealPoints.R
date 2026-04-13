# --- Input validation tests ---

test_that("mapBackOnRealPoints rejects NULL dataset", {
  point <- data.frame(x = 0, y = 0)
  expect_error(mapBackOnRealPoints(dataset = NULL, point = point, dim = "x"),
               "'dataset' must be provided")
})

test_that("mapBackOnRealPoints rejects NULL point", {
  dataset <- data.frame(x = 1:5, y = 1:5)
  expect_error(mapBackOnRealPoints(dataset = dataset, point = NULL, dim = "x"),
               "'point' must be provided")
})

test_that("mapBackOnRealPoints rejects non-dataframe dataset", {
  expect_error(mapBackOnRealPoints(dataset = "bad", point = data.frame(x = 1), dim = "x"),
               "'dataset' must be a data.frame")
})

test_that("mapBackOnRealPoints rejects empty dataset", {
  expect_error(
    mapBackOnRealPoints(dataset = data.frame(x = numeric(0)), point = data.frame(x = 1), dim = "x"),
    "at least one row"
  )
})

test_that("mapBackOnRealPoints rejects bad dim", {
  dataset <- data.frame(x = 1:5, y = 1:5)
  point <- data.frame(x = 0, y = 0)
  expect_error(mapBackOnRealPoints(dataset = dataset, point = point, dim = ""),
               "'dim' must be a character vector")
})

test_that("mapBackOnRealPoints rejects dim not in dataset columns", {
  dataset <- data.frame(x = 1:5, y = 1:5)
  point <- data.frame(x = 0, y = 0, z = 0)
  expect_error(mapBackOnRealPoints(dataset = dataset, point = point, dim = c("x", "z")),
               "columns not found in 'dataset'")
})

test_that("mapBackOnRealPoints rejects dim not in point columns", {
  dataset <- data.frame(x = 1:5, y = 1:5, z = 1:5)
  point <- data.frame(x = 0, y = 0)
  expect_error(mapBackOnRealPoints(dataset = dataset, point = point, dim = c("x", "z")),
               "columns not found in 'point'")
})

test_that("mapBackOnRealPoints rejects negative threshold", {
  dataset <- data.frame(x = 1:5, y = 1:5)
  point <- data.frame(x = 0, y = 0)
  expect_error(mapBackOnRealPoints(dataset = dataset, point = point, dim = "x", threshold = -1),
               "'threshold' must be a positive number")
})

# --- Positive tests ---

test_that("mapBackOnRealPoints finds exact match", {
  dataset <- data.frame(x = c(1, 2, 3), y = c(10, 20, 30))
  point <- data.frame(x = 2, y = 20)
  result <- mapBackOnRealPoints(dataset = dataset, point = point, dim = c("x", "y"), threshold = 1)
  expect_false(is.na(result[1]))
  expect_equal(result$distanceFromSampledPoint, 0)
})

test_that("mapBackOnRealPoints returns NA when point is too far", {
  dataset <- data.frame(x = c(100, 200), y = c(100, 200))
  point <- data.frame(x = 0, y = 0)
  result <- mapBackOnRealPoints(dataset = dataset, point = point, dim = c("x", "y"), threshold = 1)
  expect_true(is.na(result))
})

test_that("mapBackOnRealPoints returns closest point with distance", {
  dataset <- data.frame(x = c(0, 10, 20), y = c(0, 0, 0))
  point <- data.frame(x = 9, y = 0)
  result <- mapBackOnRealPoints(dataset = dataset, point = point, dim = c("x", "y"), threshold = 5)
  # Closest is (10, 0), Manhattan distance = |9-10| + |0-0| = 1
  expect_equal(result$x, 10)
  expect_equal(result$distanceFromSampledPoint, 1)
})
