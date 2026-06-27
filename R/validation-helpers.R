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

# Registry of method-specific pseudo-absence-sampler parameters, used by
# check_cross_sampler_args() to turn "argument passed to the wrong sampler" into
# actionable guidance. Keyed by parameter name; each entry has $owners (the samplers
# for which the parameter is valid) and $analog (named by each OTHER sampler -> the
# one-line instruction shown when the parameter reaches that sampler). The shared
# arguments env.rast/pres/plot_proc/verbose are real formals of every sampler and so
# never reach `...`; they are intentionally absent here.
.sampler_arg_registry <- local({
  reg <- list()
  add <- function(params, owners, analog) {
    for (p in params) reg[[p]] <<- list(owners = owners, analog = analog)
  }
  add(c("thres", "H", "grid.res", "n.tr", "prev"),
      owners = c("paSampling", "paSamplingNn"),
      analog = c(paSamplingMcmc = paste0(
        "paSamplingMcmc() has no sampling grid; size the output with 'n.samples'/",
        "'chain.length' and tune presence exclusion with 'species.cutoff.threshold'.")))
  add(c("dimensions", "precomputed.pca", "n.samples"),
      owners = c("paSamplingNn", "paSamplingMcmc"),
      analog = c(paSampling = paste0(
        "paSampling() is 2-D only and sizes output per grid cell; use 'grid.res'/'n.tr'/",
        "'prev' (and switch to paSamplingNn()/paSamplingMcmc() for >2 dimensions).")))
  add(c("nn.based.presence.exclusion", "data.based.distance.threshold",
        "n.candidates", "dim.correction"),
      owners = "paSamplingNn",
      analog = c(
        paSampling = paste0(
          "paSampling() uses a KDE + convex-hull filter on a 2-D grid; tune it with ",
          "'thres'/'H'/'grid.res'."),
        paSamplingMcmc = paste0(
          "paSamplingMcmc() samples with a Markov chain; tune exclusion with ",
          "'species.cutoff.threshold' and coverage with 'chain.length'/'n.samples'.")))
  add(c("chain.length", "burnIn", "covariance.correction", "seed.number", "num.chains",
        "num.cores", "engine", "precomputed.env", "n.neighbors.for.statistics",
        "low.end.of.inclueded.points", "high.end.of.included.points",
        "environmental.cutof.percentile"),
      owners = "paSamplingMcmc",
      analog = c(
        paSampling = paste0(
          "paSampling() has no Markov chain; control sampling with 'thres'/'grid.res'/'n.tr'."),
        paSamplingNn = paste0(
          "paSamplingNn() has no Markov chain; control coverage with 'n.candidates'/'n.tr' ",
          "and exclusion with 'thres'/'dim.correction'.")))
  add("species.cutoff.threshold",
      owners = "paSamplingMcmc",
      analog = c(
        paSampling   = "the paSampling() analog is 'thres' (KDE presence-exclusion quantile).",
        paSamplingNn = "the paSamplingNn() analog is 'thres' (presence-exclusion quantile)."))
  reg
})

#' Reject (or guide) arguments passed to the wrong pseudo-absence sampler
#'
#' Each sampler captures unmatched arguments via \code{...} and calls this helper as its
#' first step, so that passing another sampler's parameter yields actionable guidance
#' instead of R's opaque "unused argument" error. The deprecated \code{env.data.raster}
#' name maps to \code{env.rast}; a name belonging to a different sampler reports the
#' analogous parameter; a name in neither the registry nor the formals (a typo) is
#' rejected outright. Errors use \code{stop_config()} so wrappers that \code{tryCatch}
#' the samplers re-raise them loudly instead of silently falling back.
#' @param dots List of the \code{...} arguments captured by the sampler.
#' @param this.sampler character(1): the sampler that was called.
#' @keywords internal
check_cross_sampler_args <- function(dots, this.sampler) {
  nms <- names(dots)
  if (is.null(nms) || !length(nms)) return(invisible(NULL))
  for (nm in nms) {
    if (!nzchar(nm)) {
      stop_config("unexpected unnamed argument passed to ", this.sampler, "().")
    }
    if (identical(nm, "env.data.raster")) {
      stop_config("'env.data.raster' has been renamed to 'env.rast'. Pass env.rast = ... instead.")
    }
    entry <- .sampler_arg_registry[[nm]]
    if (is.null(entry)) {
      stop_config("unknown argument '", nm, "' passed to ", this.sampler, "().")
    }
    if (!this.sampler %in% entry$owners) {
      stop_config("'", nm, "' is a ", paste(entry$owners, collapse = "/"),
                  "() parameter, not ", this.sampler, "(). ", entry$analog[[this.sampler]])
    }
  }
  invisible(NULL)
}
