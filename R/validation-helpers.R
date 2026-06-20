# Internal validation helper functions for USE.MCMC
# These are not exported and are used across multiple functions
# to provide consistent, informative error messages.

#' Check that a raster input is valid
#' @param x The object to check
#' @param arg_name Name of the argument (for error messages)
#' @keywords internal
check_raster_input <- function(x, arg_name = "env.rast") {
  if (is.null(x)) {
    stop(paste0("'", arg_name, "' must be provided (got NULL)"), call. = FALSE)
  }
  if (!inherits(x, "BasicRaster") && !inherits(x, "SpatRaster")) {
    stop(paste0("'", arg_name, "' must be a SpatRaster or BasicRaster object, got '",
                paste(class(x), collapse = "/"), "'"), call. = FALSE)
  }
}

#' Check that a spatial points input is valid
#' @param x The object to check
#' @param arg_name Name of the argument (for error messages)
#' @param allow_null If TRUE, NULL values are accepted without error
#' @keywords internal
check_spatial_points <- function(x, arg_name = "pres", allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(invisible(NULL))
    stop(paste0("'", arg_name, "' must be provided (got NULL)"), call. = FALSE)
  }
  if (!inherits(x, "SpatialPoints") && !inherits(x, "SpatialPointsDataFrame") &&
      !inherits(x, "SpatVector") && !inherits(x, "sf")) {
    stop(paste0("'", arg_name, "' must be a spatial object (sf, SpatVector, SpatialPoints, ",
                "or SpatialPointsDataFrame), got '", paste(class(x), collapse = "/"), "'"),
         call. = FALSE)
  }
}

#' Check that column names exist in a dataset
#' @param data The data.frame or sf object to check columns in
#' @param columns Character vector of column names to look for
#' @param data_name Name of the data argument (for error messages)
#' @param col_name Name of the columns argument (for error messages)
#' @keywords internal
check_columns_exist <- function(data, columns, data_name = "dataset", col_name = "dimensions") {
  if (inherits(data, "sf")) {
    available <- setdiff(names(data), attr(data, "sf_column"))
  } else {
    available <- names(data)
  }
  missing_cols <- setdiff(columns, available)
  if (length(missing_cols) > 0) {
    stop(paste0("'", col_name, "' contains columns not found in '", data_name, "': ",
                paste0("'", missing_cols, "'", collapse = ", "),
                ". Available columns: ", paste(available, collapse = ", ")),
         call. = FALSE)
  }
}

#' Check that a value is a positive integer
#' @param x The value to check
#' @param arg_name Name of the argument (for error messages)
#' @param allow_zero If TRUE, zero is accepted
#' @keywords internal
check_positive_integer <- function(x, arg_name, allow_zero = FALSE) {
  if (is.null(x)) {
    stop(paste0("'", arg_name, "' must be provided (got NULL)"), call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1) {
    stop(paste0("'", arg_name, "' must be a single numeric value, got '",
                paste(class(x), collapse = "/"), "' of length ", length(x)),
         call. = FALSE)
  }
  if (is.na(x)) {
    stop(paste0("'", arg_name, "' must not be NA"), call. = FALSE)
  }
  if (!allow_zero && x <= 0) {
    stop(paste0("'", arg_name, "' must be a positive integer, got ", x), call. = FALSE)
  }
  if (allow_zero && x < 0) {
    stop(paste0("'", arg_name, "' must be a non-negative integer, got ", x), call. = FALSE)
  }
  if (x != floor(x)) {
    stop(paste0("'", arg_name, "' must be an integer value, got ", x), call. = FALSE)
  }
}

#' Check that numeric values fall within a range
#' @param x The value(s) to check
#' @param arg_name Name of the argument (for error messages)
#' @param min_val Minimum allowed value (inclusive)
#' @param max_val Maximum allowed value (inclusive)
#' @keywords internal
check_in_range <- function(x, arg_name, min_val = 0, max_val = 1) {
  if (is.null(x)) {
    stop(paste0("'", arg_name, "' must be provided (got NULL)"), call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(paste0("'", arg_name, "' must be numeric, got '",
                paste(class(x), collapse = "/"), "'"), call. = FALSE)
  }
  if (any(is.na(x))) {
    stop(paste0("'", arg_name, "' must not contain NA values"), call. = FALSE)
  }
  out_of_range <- x[x < min_val | x > max_val]
  if (length(out_of_range) > 0) {
    stop(paste0("All '", arg_name, "' values must be between ", min_val, " and ", max_val,
                ". Got values: ", paste(out_of_range, collapse = ", ")),
         call. = FALSE)
  }
}

#' Signal a configuration error (bad arguments / incompatible cache)
#'
#' Raises a condition of class \code{"USE.MCMC_config_error"} (in addition to the
#' usual \code{"error"}). Callers that wrap paSamplingMcmc() in a tryCatch and
#' fall back to another sampler on failure can detect this class and re-raise it,
#' so that a misconfiguration (e.g. an incompatible \code{precomputed.env}) fails
#' loudly instead of being silently downgraded to a different sampler.
#' @param ... Parts of the error message, pasted together.
#' @keywords internal
stop_config <- function(...) {
  stop(structure(
    class = c("USE.MCMC_config_error", "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}
