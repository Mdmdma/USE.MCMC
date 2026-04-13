test_that("alwaysOne returns 1 with no arguments", {
  expect_equal(alwaysOne(), 1)
})

test_that("alwaysOne returns 1 with multiple arguments", {
  expect_equal(alwaysOne(1, 2, 3), 1)
})

test_that("alwaysOne returns 1 with NULL input", {
  expect_equal(alwaysOne(NULL), 1)
})

test_that("alwaysOne returns 1 even on error input", {
  expect_equal(alwaysOne(stop("this should be an error")), 1)
})

test_that("alwaysOne returns numeric of length 1", {
  result <- alwaysOne()
  expect_true(is.numeric(result))
  expect_length(result, 1)
})
