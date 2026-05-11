test_that("acceptNextPoint accepts when proposed density is higher", {
  expect_true(USE.MCMC:::acceptNextPoint(1, 5))
})

test_that("acceptNextPoint accepts when densities are equal", {
  # ratio = 1, so acceptance.ratio > 1 is FALSE
  # but runif(1) < 1 is almost surely TRUE
  set.seed(42)
  expect_true(USE.MCMC:::acceptNextPoint(1, 1))
})

test_that("acceptNextPoint mostly rejects much lower density", {
  # ratio = 0.00001, very unlikely to accept
  set.seed(42)
  results <- replicate(50, USE.MCMC:::acceptNextPoint(100, 0.001))
  # Should reject most of the time
  expect_true(sum(results) < 10)
})

test_that("acceptNextPoint rejects when current density is zero (degenerate)", {
  # Non-positive current density signals an unreachable starting point;
  # any proposal is rejected to keep R and C++ engines in lockstep.
  expect_false(USE.MCMC:::acceptNextPoint(0, 1))
})

test_that("acceptNextPoint rejects when proposed density is zero or negative", {
  expect_false(USE.MCMC:::acceptNextPoint(1, 0))
  expect_false(USE.MCMC:::acceptNextPoint(1, -0.5))
})

test_that("acceptNextPoint rejects when both densities are zero", {
  expect_false(USE.MCMC:::acceptNextPoint(0, 0))
})

test_that("acceptNextPoint rejects on NA densities", {
  expect_false(USE.MCMC:::acceptNextPoint(NA_real_, 1))
  expect_false(USE.MCMC:::acceptNextPoint(1, NA_real_))
})
