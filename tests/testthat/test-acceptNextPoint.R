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

test_that("acceptNextPoint accepts when current density is zero", {
  # ratio = 1/0 = Inf, Inf > 1 is TRUE
  expect_true(USE.MCMC:::acceptNextPoint(0, 1))
})

test_that("acceptNextPoint rejects when both densities are zero (NaN ratio)", {
  # ratio = 0/0 = NaN, handled gracefully by rejecting the proposal
  expect_false(USE.MCMC:::acceptNextPoint(0, 0))
})
