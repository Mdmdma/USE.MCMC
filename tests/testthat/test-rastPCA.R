# --- Input validation tests ---

test_that("rastPCA rejects non-raster input", {
  expect_error(rastPCA("not a raster"), "must be a SpatRaster or BasicRaster")
})

test_that("rastPCA rejects NULL input", {
  expect_error(rastPCA(NULL), "must be provided")
})

test_that("rastPCA rejects single-layer raster", {
  r <- make_test_raster()
  expect_error(rastPCA(r[[1]]), "At least two layers")
})

test_that("rastPCA rejects non-logical naMask", {
  r <- make_test_raster()
  expect_error(rastPCA(r, naMask = "yes"), "'naMask' must be a single logical")
})

test_that("rastPCA rejects non-logical stand", {
  r <- make_test_raster()
  expect_error(rastPCA(r, stand = 1), "'stand' must be a single logical")
})

test_that("rastPCA rejects non-integer nPC", {
  r <- make_test_raster()
  expect_error(rastPCA(r, nPC = 1.5), "'nPC' must be a positive integer")
})

test_that("rastPCA rejects negative nPC", {
  r <- make_test_raster()
  expect_error(rastPCA(r, nPC = -1), "'nPC' must be a positive integer")
})

test_that("rastPCA messages when nPC exceeds layers", {
  r <- make_test_raster()
  expect_message(rastPCA(r, nPC = 100), "maximum number of PCs")
})

# --- Positive tests ---

test_that("rastPCA returns correct structure", {
  r <- make_test_raster()
  result <- rastPCA(r)
  expect_type(result, "list")
  expect_true("call" %in% names(result))
  expect_true("pca" %in% names(result))
  expect_true("PCs" %in% names(result))
})

test_that("rastPCA pca element is princomp", {
  r <- make_test_raster()
  result <- rastPCA(r)
  expect_s3_class(result$pca, "princomp")
})

test_that("rastPCA PCs element is SpatRaster", {
  r <- make_test_raster()
  result <- rastPCA(r)
  expect_true(inherits(result$PCs, "SpatRaster"))
})

test_that("rastPCA returns correct number of PCs", {
  r <- make_test_raster()
  nlayers <- terra::nlyr(r)
  result <- rastPCA(r)
  expect_equal(terra::nlyr(result$PCs), nlayers)
})

test_that("rastPCA returns requested number of PCs", {
  r <- make_test_raster()
  result <- rastPCA(r, nPC = 2)
  expect_equal(terra::nlyr(result$PCs), 2)
})

test_that("rastPCA PC names follow pattern", {
  r <- make_test_raster()
  result <- rastPCA(r, nPC = 3)
  expect_equal(names(result$PCs), c("PC1", "PC2", "PC3"))
})

test_that("rastPCA works with stand=TRUE", {
  r <- make_test_raster()
  result <- rastPCA(r, stand = TRUE)
  expect_s3_class(result$pca, "princomp")
  expect_true(!is.null(result$pca$scale))
})
