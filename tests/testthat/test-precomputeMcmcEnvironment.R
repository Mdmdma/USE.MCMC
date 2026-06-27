skip_on_cran()

# Shared fixtures (file-level code runs once). Small chain so the suite is quick
# while still exercising the real MCMC path.
r        <- make_test_raster()
pres     <- make_test_presence_sf(r, 30)
dims     <- c("PC1", "PC2")
seed     <- 7L
run_args <- list(n.samples = 50, chain.length = 2000, burnIn = 200,
                 num.chains = 1, num.cores = 1, dimensions = dims, seed.number = seed)

env_bundle <- precomputeMcmcEnvironment(env.rast = r, dimensions = dims, seed.number = seed)

do_run <- function(...) {
  do.call(paSamplingMcmc, c(list(env.rast = r, pres = pres), run_args, list(...)))
}
ref_uncached <- do_run()

test_that("precomputeMcmcEnvironment returns a complete, serialisable bundle", {
  expect_s3_class(env_bundle, "mcmc_environment")
  expect_true(all(c("pcs_packed", "env.with.pc.sf", "env.data.cleaned",
                    "environmental.data.model", "environmental.densities",
                    "covariance.matrix", "distance.threshold", "dimensions",
                    "rng_state") %in% names(env_bundle)))
  expect_identical(env_bundle$dimensions, dims)
  expect_s4_class(terra::unwrap(env_bundle$pcs_packed), "SpatRaster")
})

test_that("precomputeMcmcEnvironment validates its inputs", {
  expect_error(precomputeMcmcEnvironment(env.rast = NULL), "must be provided")
  expect_error(precomputeMcmcEnvironment(env.rast = r, dimensions = "PC1"),
               "at least 2 elements")
  expect_error(precomputeMcmcEnvironment(env.rast = r, seed.number = "x"),
               "seed.number")
})

test_that("cached result is bit-identical to the uncached path", {
  res_cached <- do_run(precomputed.env = env_bundle)
  expect_equal(sf::st_coordinates(res_cached), sf::st_coordinates(ref_uncached))
})

test_that("the bundle survives a saveRDS/readRDS round-trip", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(env_bundle, f)
  res_rt <- do_run(precomputed.env = readRDS(f))
  expect_equal(sf::st_coordinates(res_rt), sf::st_coordinates(ref_uncached))
})

test_that("env.rast is optional once a bundle is supplied", {
  res <- do.call(paSamplingMcmc,
                 c(list(env.rast = NULL, pres = pres), run_args,
                   list(precomputed.env = env_bundle)))
  expect_equal(sf::st_coordinates(res), sf::st_coordinates(ref_uncached))
})

test_that("an incompatible precomputed.env raises a config error (not a silent fallback)", {
  bad_missing <- env_bundle; bad_missing$rng_state <- NULL
  expect_error(do_run(precomputed.env = bad_missing), class = "USE.MCMC_config_error")
  expect_error(do_run(precomputed.env = list(a = 1)), class = "USE.MCMC_config_error")
})

test_that("a dimension mismatch raises a config error", {
  wrong_dims <- env_bundle; wrong_dims$dimensions <- c("PC1", "PC2", "PC3")
  expect_error(do_run(precomputed.env = wrong_dims), class = "USE.MCMC_config_error")
})
