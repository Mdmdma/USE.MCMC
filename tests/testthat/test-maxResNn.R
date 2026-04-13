skip_on_cran()

# --- Input validation tests ---

test_that("maxResNn rejects non-raster input when PCA=FALSE", {
  expect_error(maxResNn("bad"), "must be a SpatRaster or BasicRaster")
})

test_that("maxResNn rejects non-logical PCA", {
  expect_error(maxResNn(NULL, PCA = "yes"), "'PCA' must be a single logical")
})

test_that("maxResNn rejects invalid PCA=TRUE input", {
  expect_error(maxResNn(list(a = 1), PCA = TRUE), "must be a list with a '\\$PCs' element")
})

test_that("maxResNn rejects non-character dimensions", {
  r <- make_test_raster()
  expect_error(maxResNn(r, dimensions = c(1, 2)), "'dimensions' must be a character vector")
})

test_that("maxResNn rejects single dimension", {
  r <- make_test_raster()
  expect_error(maxResNn(r, dimensions = "PC1"), "at least 2 elements")
})

test_that("maxResNn rejects non-positive n.neighbors", {
  r <- make_test_raster()
  expect_error(maxResNn(r, n.neighbors = 0), "'n.neighbors' must be a positive integer")
})

# --- Positive tests ---

test_that("maxResNn returns positive numeric", {
  r <- make_test_raster()
  result <- maxResNn(r, dimensions = c("PC1", "PC2"))
  expect_true(is.numeric(result))
  expect_true(result > 0)
})

test_that("maxResNn works with PCA=TRUE and precomputed PCA", {
  r <- make_test_raster()
  rpc <- rastPCA(r, stand = TRUE)
  result <- maxResNn(rpc, dimensions = c("PC1", "PC2"), PCA = TRUE)
  expect_true(is.numeric(result))
  expect_true(result > 0)
})
