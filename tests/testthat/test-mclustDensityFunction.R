# --- Input validation tests ---

test_that("mclustDensityFunction rejects NULL env.model", {
  expect_error(mclustDensityFunction(env.model = NULL),
               "'env.model' must be provided")
})

test_that("mclustDensityFunction rejects non-densityMclust env.model", {
  expect_error(mclustDensityFunction(env.model = list(a = 1), species.model = NULL),
               "'env.model' must be a densityMclust object")
})

test_that("mclustDensityFunction accepts NULL species.model (uniform mode)", {
  # species.model = NULL is the uniform-environment target: 1 inside the support,
  # floor outside, no presence GMM subtracted.
  set.seed(42)
  env_data <- data.frame(x = rnorm(50), y = rnorm(50))
  env_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  fn <- mclustDensityFunction(env.model = env_model, species.model = NULL,
                              dim = c("x", "y"), threshold = 0.01)
  expect_true(is.function(fn))
  # Inside the environmental support the uniform target is exactly 1.
  expect_identical(fn(c(0, 0)), 1)
  # Far outside the support it drops to the floor (threshold / 1000).
  expect_true(fn(c(100, 100)) < 0.01)
  # The Rcpp spec carries the infinite-cutoff sentinel so the C++ inner loop
  # returns the same uniform target without a dedicated branch.
  spec <- attr(fn, "rcpp_spec")
  expect_identical(spec$type, "mclust_density")
  expect_true(is.infinite(spec$species_cutoff))
})

test_that("mclustDensityFunction rejects non-densityMclust species.model", {
  set.seed(42)
  env_data <- data.frame(x = rnorm(50), y = rnorm(50))
  env_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  expect_error(mclustDensityFunction(env.model = env_model, species.model = "bad"),
               "'species.model' must be a densityMclust object")
})

test_that("mclustDensityFunction rejects bad dim", {
  set.seed(42)
  env_data <- data.frame(x = rnorm(50), y = rnorm(50))
  env_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  sp_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  expect_error(mclustDensityFunction(env.model = env_model, species.model = sp_model, dim = 123),
               "'dim' must be a character vector")
})

test_that("mclustDensityFunction rejects non-positive threshold", {
  set.seed(42)
  env_data <- data.frame(x = rnorm(50), y = rnorm(50))
  env_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  sp_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  expect_error(mclustDensityFunction(env.model = env_model, species.model = sp_model,
                                     dim = c("x", "y"), threshold = -1),
               "'threshold' must be a positive number")
})

test_that("mclustDensityFunction rejects non-positive species.cutoff.threshold", {
  set.seed(42)
  env_data <- data.frame(x = rnorm(50), y = rnorm(50))
  env_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  sp_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  expect_error(mclustDensityFunction(env.model = env_model, species.model = sp_model,
                                     dim = c("x", "y"), species.cutoff.threshold = 0),
               "'species.cutoff.threshold' must be a positive number")
})

# --- Positive tests ---

test_that("mclustDensityFunction returns a function", {
  set.seed(42)
  env_data <- data.frame(x = rnorm(50), y = rnorm(50))
  env_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  sp_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  fn <- mclustDensityFunction(env.model = env_model, species.model = sp_model,
                              dim = c("x", "y"))
  expect_true(is.function(fn))
})

test_that("returned density function returns numeric", {
  set.seed(42)
  env_data <- data.frame(x = rnorm(50), y = rnorm(50))
  env_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  sp_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  fn <- mclustDensityFunction(env.model = env_model, species.model = sp_model,
                              dim = c("x", "y"), threshold = 0.01)
  point <- c(0, 0)
  result <- fn(point)
  expect_true(is.numeric(result))
  expect_length(result, 1)
})

test_that("density function returns low value for far-out points", {
  set.seed(42)
  env_data <- data.frame(x = rnorm(100), y = rnorm(100))
  env_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  sp_model <- mclust::densityMclust(env_data, plot = FALSE, verbose = FALSE)
  fn <- mclustDensityFunction(env.model = env_model, species.model = sp_model,
                              dim = c("x", "y"), threshold = 0.01)
  # Point very far from the data center, should be below threshold
  far_point <- c(100, 100)
  result <- fn(far_point)
  expect_true(result <= 0.01)
})
