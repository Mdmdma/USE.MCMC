#' Map Back on Real Points
#' searches the closest point in the dataset regarding the given point and the dimensions given
#'
#' @param dataset dataframe of the target dataset
#' @param point dataframe containing the given point
#' @param dim string vector containing the names of the dimensions that should be respected
#' @param threshold threshold obove which we consider the distance too big and want to discard the point
#'
#' @returns closest point that we could find in the dataset with the distance added as a parameter
#' @export
#'
mapBackOnRealPoints <- function(dataset = NULL, point = NULL, dim ="", threshold = 0.5){
  # Input validation
  if (is.null(dataset)) {
    stop("'dataset' must be provided (got NULL)", call. = FALSE)
  }
  if (!is.data.frame(dataset) && !inherits(dataset, "sf")) {
    stop(paste0("'dataset' must be a data.frame or sf object, got '",
                paste(class(dataset), collapse = "/"), "'"), call. = FALSE)
  }
  if (nrow(dataset) == 0) {
    stop("'dataset' must have at least one row", call. = FALSE)
  }
  if (is.null(point)) {
    stop("'point' must be provided (got NULL)", call. = FALSE)
  }
  if (!is.data.frame(point) && !inherits(point, "sf")) {
    stop(paste0("'point' must be a data.frame or sf object, got '",
                paste(class(point), collapse = "/"), "'"), call. = FALSE)
  }
  if (!is.character(dim) || length(dim) < 1 || all(dim == "")) {
    stop("'dim' must be a character vector of column names", call. = FALSE)
  }
  # Check columns exist in both dataset and point
  ds_cols <- if (inherits(dataset, "sf")) setdiff(names(dataset), attr(dataset, "sf_column")) else names(dataset)
  pt_cols <- if (inherits(point, "sf")) setdiff(names(point), attr(point, "sf_column")) else names(point)
  missing_ds <- setdiff(dim, ds_cols)
  missing_pt <- setdiff(dim, pt_cols)
  if (length(missing_ds) > 0) {
    stop(paste0("'dim' contains columns not found in 'dataset': ",
                paste0("'", missing_ds, "'", collapse = ", "),
                ". Available: ", paste(ds_cols, collapse = ", ")), call. = FALSE)
  }
  if (length(missing_pt) > 0) {
    stop(paste0("'dim' contains columns not found in 'point': ",
                paste0("'", missing_pt, "'", collapse = ", "),
                ". Available: ", paste(pt_cols, collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0) {
    stop(paste0("'threshold' must be a positive number, got ", deparse(threshold)), call. = FALSE)
  }

  distances <- apply(dataset, 1, function(row) euclidianMetric(row, point, dim))
  min.dist.index <- which.min(distances)
  min.dist <- distances[min.dist.index]
  best.point <- dataset[min.dist.index,]
  best.point$distanceFromSampledPoint <- min.dist
  if (min.dist < threshold) return(best.point)
  return(NA)
}
