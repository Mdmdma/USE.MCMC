skip_on_cran()

# --- Input validation tests ---

test_that("mcmcSampling rejects NULL dataset", {
  expect_error(mcmcSampling(dataset = NULL), "'dataset' must be provided")
})

test_that("mcmcSampling rejects non-dataframe dataset", {
  expect_error(mcmcSampling(dataset = "bad"), "'dataset' must be a data.frame")
})

test_that("mcmcSampling rejects empty dataset", {
  expect_error(
    mcmcSampling(dataset = data.frame(x = numeric(0))),
    "at least one row"
  )
})

test_that("mcmcSampling rejects dimensions not in dataset", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    mcmcSampling(dataset = df, dimensions = c("x", "z")),
    "columns not found in 'dataset'"
  )
})

test_that("mcmcSampling rejects non-function densityFunction", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    mcmcSampling(dataset = df, dimensions = c("x", "y"),
                 densityFunction = "not a function",
                 proposalFunction = addHighDimGaussian(dim = 2)),
    "'densityFunction' must be a function"
  )
})

test_that("mcmcSampling rejects non-function proposalFunction", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    mcmcSampling(dataset = df, dimensions = c("x", "y"),
                 proposalFunction = 42),
    "'proposalFunction' must be a function"
  )
})

test_that("mcmcSampling rejects negative n.sample.points", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    mcmcSampling(dataset = df, dimensions = c("x", "y"),
                 proposalFunction = addHighDimGaussian(dim = 2),
                 n.sample.points = -1),
    "'n.sample.points' must be a positive number"
  )
})

test_that("mcmcSampling rejects negative burnIn", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    mcmcSampling(dataset = df, dimensions = c("x", "y"),
                 proposalFunction = addHighDimGaussian(dim = 2),
                 n.sample.points = 10, burnIn = -5),
    "'burnIn' must be a non-negative number"
  )
})

test_that("mcmcSampling rejects non-logical verbose", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    mcmcSampling(dataset = df, dimensions = c("x", "y"),
                 proposalFunction = addHighDimGaussian(dim = 2),
                 n.sample.points = 10, verbose = "yes"),
    "'verbose' must be a single logical"
  )
})

test_that("mcmcSampling rejects non-positive covariance.correction", {
  df <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    mcmcSampling(dataset = df, dimensions = c("x", "y"),
                 proposalFunction = addHighDimGaussian(dim = 2),
                 n.sample.points = 10, covariance.correction = 0),
    "'covariance.correction' must be a positive number"
  )
})

# --- Positive tests ---

test_that("mcmcSampling returns dataframe with correct dimensions", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(50)
  result <- mcmcSampling(
    dataset = pc_sf,
    dimensions = c("PC1", "PC2"),
    densityFunction = alwaysOne,
    proposalFunction = addHighDimGaussian(dim = 2),
    n.sample.points = 10,
    burnIn = 0,
    verbose = FALSE,
    covariance.correction = 1
  )
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 10)
})

test_that("mcmcSampling output has dimension columns", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(50)
  result <- mcmcSampling(
    dataset = pc_sf,
    dimensions = c("PC1", "PC2"),
    densityFunction = alwaysOne,
    proposalFunction = addHighDimGaussian(dim = 2),
    n.sample.points = 10,
    burnIn = 0,
    verbose = FALSE,
    covariance.correction = 1
  )
  expect_true("PC1" %in% names(result))
  expect_true("PC2" %in% names(result))
  expect_true("density" %in% names(result))
})

test_that("mcmcSampling output has non-NA dimension values", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(50)
  result <- mcmcSampling(
    dataset = pc_sf,
    dimensions = c("PC1", "PC2"),
    densityFunction = alwaysOne,
    proposalFunction = addHighDimGaussian(dim = 2),
    n.sample.points = 10,
    burnIn = 0,
    verbose = FALSE,
    covariance.correction = 1
  )
  expect_true(all(!is.na(result$PC1)))
  expect_true(all(!is.na(result$PC2)))
})

test_that("mcmcSampling works with sf input", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(50)
  expect_no_error(
    mcmcSampling(
      dataset = pc_sf,
      dimensions = c("PC1", "PC2"),
      densityFunction = alwaysOne,
      proposalFunction = addHighDimGaussian(dim = 2),
      n.sample.points = 5,
      burnIn = 0,
      verbose = FALSE
    )
  )
})
