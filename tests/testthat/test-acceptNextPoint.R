test_that("acceptNextPoint accepts when proposed density is higher", {
  current <- data.frame(x = 0, density = 1)
  proposed <- data.frame(x = 1, density = 5)
  expect_true(USE.MCMC:::acceptNextPoint(current, proposed))
})

test_that("acceptNextPoint accepts when densities are equal", {
  current <- data.frame(x = 0, density = 1)
  proposed <- data.frame(x = 1, density = 1)
  # ratio = 1, so acceptance.ratio > 1 is FALSE

  # but runif(1) < 1 is almost surely TRUE
  set.seed(42)
  expect_true(USE.MCMC:::acceptNextPoint(current, proposed))
})

test_that("acceptNextPoint mostly rejects much lower density", {
  current <- data.frame(x = 0, density = 100)
  proposed <- data.frame(x = 1, density = 0.001)
  # ratio = 0.00001, very unlikely to accept
  set.seed(42)
  results <- replicate(50, USE.MCMC:::acceptNextPoint(current, proposed))
  # Should reject most of the time

  expect_true(sum(results) < 10)
})

test_that("acceptNextPoint accepts when current density is zero", {
  current <- data.frame(x = 0, density = 0)
  proposed <- data.frame(x = 1, density = 1)
  # ratio = 1/0 = Inf, Inf > 1 is TRUE
  expect_true(USE.MCMC:::acceptNextPoint(current, proposed))
})

test_that("acceptNextPoint rejects when both densities are zero (NaN ratio)", {
  current <- data.frame(x = 0, density = 0)
  proposed <- data.frame(x = 1, density = 0)
  # ratio = 0/0 = NaN, handled gracefully by rejecting the proposal
  expect_false(USE.MCMC:::acceptNextPoint(current, proposed))
})
