test_that("function returns 1 even on error input", {
  expect_equal(alwaysOne(stop("this should be an error")), 1)
})
