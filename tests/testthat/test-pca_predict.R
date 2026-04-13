# pca_predict is an internal function called indirectly by terra::predict
# It cannot be tested in isolation because predict(data, model) dispatches
# on the class of data. It is exercised through rastPCA which uses
# terra::predict(env.rast, eigenDecomp, nPC=nPC, fun=pca_predict).

test_that("pca_predict exists and is a function", {
  expect_true(is.function(USE.MCMC:::pca_predict))
})

test_that("pca_predict has correct formals", {
  fmls <- formals(USE.MCMC:::pca_predict)
  expect_true("data" %in% names(fmls))
  expect_true("model" %in% names(fmls))
  expect_true("nPC" %in% names(fmls))
})

test_that("pca_predict is covered via rastPCA integration", {
  # rastPCA internally calls terra::predict(..., fun = pca_predict)
  r <- make_test_raster()
  result <- rastPCA(r, nPC = 2)
  # If pca_predict works, PCs will have the right structure
  expect_equal(terra::nlyr(result$PCs), 2)
  expect_equal(names(result$PCs), c("PC1", "PC2"))
})
