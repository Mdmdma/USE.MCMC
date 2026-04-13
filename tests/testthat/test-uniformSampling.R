# --- Input validation tests ---

test_that("uniformSampling rejects non-sf input", {
  expect_error(uniformSampling(sdf = data.frame(x = 1), grid.res = 5), "'sdf' must be an sf object")
})

test_that("uniformSampling rejects sf without POINT geometry", {
  poly <- sf::st_sfc(sf::st_polygon(list(cbind(c(0,1,1,0,0), c(0,0,1,1,0)))))
  sf_poly <- sf::st_sf(geometry = poly)
  expect_error(uniformSampling(sdf = sf_poly, grid.res = 5), "POINT geometry")
})

test_that("uniformSampling rejects non-numeric n.tr", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(uniformSampling(sdf = pc_sf, grid.res = 3, n.tr = "five"), "'n.tr' must be a positive number")
})

test_that("uniformSampling rejects non-logical plot_proc", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(uniformSampling(sdf = pc_sf, grid.res = 3, plot_proc = 1), "'plot_proc' must be a single logical")
})

test_that("uniformSampling rejects non-numeric grid.res", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(uniformSampling(sdf = pc_sf, grid.res = "bad"), "'grid.res' must be a positive number")
})

test_that("uniformSampling rejects non-logical sub.ts", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(uniformSampling(sdf = pc_sf, grid.res = 3, sub.ts = "yes"), "'sub.ts' must be a single logical")
})

test_that("uniformSampling rejects non-logical verbose", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(uniformSampling(sdf = pc_sf, grid.res = 3, verbose = 0), "'verbose' must be a single logical")
})

# --- Positive tests ---

test_that("uniformSampling returns sf object", {
  pc_sf <- make_test_pc_sf(200)
  result <- uniformSampling(sdf = pc_sf, grid.res = 3, n.tr = 2, verbose = FALSE, plot_proc = FALSE)
  expect_true(inherits(result, "sf"))
})

test_that("uniformSampling with sub.ts returns list with obs.tr and obs.ts", {
  pc_sf <- make_test_pc_sf(200)
  result <- uniformSampling(sdf = pc_sf, grid.res = 3, n.tr = 2, sub.ts = TRUE,
                            n.ts = 1, verbose = FALSE, plot_proc = FALSE)
  expect_type(result, "list")
  expect_true("obs.tr" %in% names(result))
  expect_true("obs.ts" %in% names(result))
  expect_true(inherits(result$obs.tr, "sf"))
  expect_true(inherits(result$obs.ts, "sf"))
})

test_that("uniformSampling result has no duplicate IDs", {
  pc_sf <- make_test_pc_sf(200)
  result <- uniformSampling(sdf = pc_sf, grid.res = 3, n.tr = 2, verbose = FALSE, plot_proc = FALSE)
  expect_false(any(duplicated(result$ID)))
})
