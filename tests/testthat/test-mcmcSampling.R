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

# --- engine = "cpp" path -----------------------------------------------------

# Helper: build a built-in density + proposal so the cpp dispatch fires.
make_cpp_fixtures <- function(pc_sf) {
  df <- sf::st_drop_geometry(pc_sf[c("PC1", "PC2")])
  env.model <- mclust::densityMclust(df, plot = FALSE, verbose = FALSE)
  sp.model <- mclust::densityMclust(df[seq_len(min(30, nrow(df))), ],
                                    plot = FALSE, verbose = FALSE)
  dfn <- mclustDensityFunction(
    env.model = env.model, species.model = sp.model,
    dim = c("PC1", "PC2"),
    threshold = stats::quantile(env.model$density, 0.001),
    species.cutoff.threshold = stats::quantile(sp.model$density, 0.95)
  )
  pfn <- addHighDimGaussian(cov.mat = stats::cov(df), dim = 2)
  list(densityFunction = dfn, proposalFunction = pfn)
}

test_that("engine = 'cpp' produces the same output shape as 'R'", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(200)
  fx <- make_cpp_fixtures(pc_sf)
  set.seed(1)
  r_out <- mcmcSampling(dataset = pc_sf, dimensions = c("PC1", "PC2"),
                        densityFunction = fx$densityFunction,
                        proposalFunction = fx$proposalFunction,
                        n.sample.points = 2000, burnIn = 500,
                        verbose = FALSE, engine = "R")
  set.seed(1)
  c_out <- mcmcSampling(dataset = pc_sf, dimensions = c("PC1", "PC2"),
                        densityFunction = fx$densityFunction,
                        proposalFunction = fx$proposalFunction,
                        n.sample.points = 2000, burnIn = 500,
                        verbose = FALSE, engine = "cpp")
  expect_equal(nrow(r_out), nrow(c_out))
  expect_equal(colnames(r_out), colnames(c_out))
  # Distributional sanity. Seeds don't yield bit-identical chains (RNG is
  # consumed in different orders), and MCMC autocorrelation makes the naive
  # sd/sqrt(N) SE under-estimate the true uncertainty in the mean. Check:
  #   marginal sd within 15% (sd is robust to autocorrelation)
  #   marginal mean within 0.3 sd of each other (loose enough to absorb
  #     correlated noise at N=500, tight enough to catch a wrong-center bug)
  expect_lt(abs(sd(r_out$PC1) - sd(c_out$PC1)) / sd(r_out$PC1), 0.15)
  expect_lt(abs(sd(r_out$PC2) - sd(c_out$PC2)) / sd(r_out$PC2), 0.15)
  expect_lt(abs(mean(r_out$PC1) - mean(c_out$PC1)) / sd(r_out$PC1), 0.3)
  expect_lt(abs(mean(r_out$PC2) - mean(c_out$PC2)) / sd(r_out$PC2), 0.3)
})

test_that("engine = 'auto' dispatches to cpp when built-in factories are used", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(200)
  fx <- make_cpp_fixtures(pc_sf)
  expect_true(!is.null(attr(fx$densityFunction, "rcpp_spec")))
  expect_true(!is.null(attr(fx$proposalFunction, "rcpp_spec")))
  # If auto dispatch wired correctly, this run completes; mirrors the cpp call.
  expect_no_error(
    mcmcSampling(dataset = pc_sf, dimensions = c("PC1", "PC2"),
                 densityFunction = fx$densityFunction,
                 proposalFunction = fx$proposalFunction,
                 n.sample.points = 50, burnIn = 10,
                 verbose = FALSE, engine = "auto")
  )
})

test_that("engine = 'cpp' errors when custom density closure has no spec", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(50)
  expect_error(
    mcmcSampling(dataset = pc_sf, dimensions = c("PC1", "PC2"),
                 densityFunction = alwaysOne,
                 proposalFunction = addHighDimGaussian(dim = 2),
                 n.sample.points = 10, burnIn = 0,
                 verbose = FALSE, engine = "cpp"),
    "engine = 'cpp' requires"
  )
})

test_that("engine = 'cpp' errors clearly when proposal covariance is singular", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(200)
  df <- sf::st_drop_geometry(pc_sf[c("PC1", "PC2")])
  env.model <- mclust::densityMclust(df, plot = FALSE, verbose = FALSE)
  sp.model <- mclust::densityMclust(df[seq_len(30), ], plot = FALSE, verbose = FALSE)
  dfn <- mclustDensityFunction(
    env.model = env.model, species.model = sp.model,
    dim = c("PC1", "PC2"),
    threshold = stats::quantile(env.model$density, 0.001),
    species.cutoff.threshold = stats::quantile(sp.model$density, 0.95)
  )
  singular_cov <- matrix(c(1, 1, 1, 1), nrow = 2)
  pfn <- addHighDimGaussian(cov.mat = singular_cov, dim = 2)
  expect_error(
    mcmcSampling(dataset = pc_sf, dimensions = c("PC1", "PC2"),
                 densityFunction = dfn, proposalFunction = pfn,
                 n.sample.points = 10, burnIn = 0,
                 verbose = FALSE, engine = "cpp"),
    "positive-definite"
  )
})

test_that("engine = 'auto' falls back to R when a custom closure is supplied", {
  set.seed(42)
  pc_sf <- make_test_pc_sf(50)
  expect_no_error(
    mcmcSampling(dataset = pc_sf, dimensions = c("PC1", "PC2"),
                 densityFunction = alwaysOne,
                 proposalFunction = addHighDimGaussian(dim = 2),
                 n.sample.points = 10, burnIn = 0,
                 verbose = FALSE, engine = "auto")
  )
})
