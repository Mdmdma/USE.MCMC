skip_on_cran()

# --- Input validation tests ---

test_that("optimRes rejects non-numeric perc.thr", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(optimRes(sdf = pc_sf, grid.res = c(2, 3), perc.thr = "ten"), "is.numeric")
})

test_that("optimRes rejects non-numeric grid.res", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(optimRes(sdf = pc_sf, grid.res = "bad"), "is.numeric")
})

test_that("optimRes rejects non-logical showOpt", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(optimRes(sdf = pc_sf, grid.res = c(2, 3), showOpt = "yes"), "is.logical")
})

test_that("optimRes rejects non-sf sdf", {
  expect_error(optimRes(sdf = data.frame(x = 1), grid.res = c(2, 3)), "sf")
})

test_that("optimRes rejects NULL cr", {
  pc_sf <- make_test_pc_sf(50)
  expect_error(optimRes(sdf = pc_sf, grid.res = c(2, 3), cr = NULL), "cluster")
})

# --- Positive tests ---

test_that("optimRes returns list with F_val and Opt_res", {
  pc_sf <- make_test_pc_sf(100)
  result <- optimRes(sdf = pc_sf, grid.res = c(2, 3, 4, 5), cr = 1, showOpt = FALSE)
  expect_type(result, "list")
  expect_true("F_val" %in% names(result))
  expect_true("Opt_res" %in% names(result))
})

test_that("optimRes F_val has correct dimensions", {
  pc_sf <- make_test_pc_sf(100)
  grid_res <- c(2, 3, 4, 5)
  result <- optimRes(sdf = pc_sf, grid.res = grid_res, cr = 1, showOpt = FALSE)
  expect_equal(nrow(result$F_val), length(grid_res))
  expect_equal(ncol(result$F_val), 2)
  expect_equal(colnames(result$F_val), c("Fval", "Res"))
})

test_that("optimRes Opt_res is numeric", {
  pc_sf <- make_test_pc_sf(100)
  result <- optimRes(sdf = pc_sf, grid.res = c(2, 3, 4, 5), cr = 1, showOpt = FALSE)
  # Opt_res can be NA if no threshold is crossed, but it should be numeric
  expect_true(is.numeric(result$Opt_res) || is.na(result$Opt_res))
})
