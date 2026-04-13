# Shared test fixtures for USE.MCMC test suite
# This file is auto-loaded by testthat before any test file runs.

# --- Minimal SpatRaster from bundled data ---
make_test_raster <- function() {
  terra::rast(USE.MCMC::Worldclim_tmp, type = "xyz")
}

# --- Minimal synthetic presence points (sf) ---
# Picks n real coordinates from the raster's non-NA cells
make_test_presence_sf <- function(env.rast, n = 20) {
  df <- terra::as.data.frame(env.rast, xy = TRUE, na.rm = TRUE)
  set.seed(42)
  idx <- sample(nrow(df), min(n, nrow(df)))
  sf::st_as_sf(df[idx, ], coords = c("x", "y"), crs = 4326)
}

# --- SpatVector variant ---
make_test_presence_vect <- function(env.rast, n = 20) {
  terra::vect(make_test_presence_sf(env.rast, n))
}

# --- Small numeric sf object for MCMC-level unit tests ---
# Two PC columns with geographic coordinates as geometry
make_test_pc_sf <- function(n = 100) {
  set.seed(42)
  df <- data.frame(
    PC1 = rnorm(n), PC2 = rnorm(n),
    x = runif(n, -10, 20), y = runif(n, 35, 55)
  )
  sf::st_as_sf(df, coords = c("x", "y"))
}

# --- Plain dataframe with PC columns ---
make_test_pc_df <- function(n = 100) {
  set.seed(42)
  data.frame(PC1 = rnorm(n), PC2 = rnorm(n))
}
