analyseAreaInHigherDim <- function(data.df, dim.to.analyse=c("PC3"), dim.to.compress=c("PC1", "PC2"), num.intervals= 20) {
  # Input validation
  if (is.null(data.df) || !is.data.frame(data.df)) {
    stop(paste0("'data.df' must be a data.frame, got '",
                paste(class(data.df), collapse = "/"), "'"), call. = FALSE)
  }
  if (nrow(data.df) == 0) {
    stop("'data.df' must have at least one row", call. = FALSE)
  }
  if (!is.character(dim.to.analyse) || length(dim.to.analyse) < 1) {
    stop("'dim.to.analyse' must be a character vector with at least 1 element", call. = FALSE)
  }
  if (!is.character(dim.to.compress) || length(dim.to.compress) < 2) {
    stop("'dim.to.compress' must be a character vector with at least 2 elements", call. = FALSE)
  }
  all_needed <- c(dim.to.analyse, dim.to.compress)
  missing_cols <- setdiff(all_needed, names(data.df))
  if (length(missing_cols) > 0) {
    stop(paste0("Columns not found in 'data.df': ",
                paste(missing_cols, collapse = ", "),
                ". Available columns: ",
                paste(names(data.df), collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(num.intervals) || length(num.intervals) != 1 || num.intervals < 1 || num.intervals != floor(num.intervals)) {
    stop(paste0("'num.intervals' must be a positive integer, got ", deparse(num.intervals)), call. = FALSE)
  }

  data.sf <- sf::st_as_sf(data.df, coords = dim.to.compress, crs = 4326)[dim.to.analyse]
  data.range <- c(min(data.df[[dim.to.analyse]], na.rm = TRUE),
                  max(data.df[[dim.to.analyse]], na.rm = TRUE))
  intervall.endpoints <- seq(from = data.range[1],
                             to = data.range[2],
                             length.out = num.intervals + 1)
  results.df <- data.frame(interval.average = numeric(num.intervals),
                           area = numeric(num.intervals))
  for (interval.counter in 1:num.intervals) {
    datapoints.in.interval <- data.sf[data.sf[[dim.to.analyse]] > intervall.endpoints[interval.counter] &
                                        data.sf[[dim.to.analyse]] < intervall.endpoints[interval.counter + 1], ]

    convex.hull <- sf::st_convex_hull(sf::st_union(datapoints.in.interval))
    area <- sf::st_area(convex.hull)
    results.df$interval.average[interval.counter]  <- mean(intervall.endpoints[interval.counter], intervall.endpoints[interval.counter + 1])
    results.df$area[interval.counter] <- area
  }
  results.df$area <- results.df$area / sum(results.df$area)
  return(results.df)
}
