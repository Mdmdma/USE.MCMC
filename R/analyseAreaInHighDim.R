analyseAreaInHigherDim <- function(data.df, dim.to.analyse=c("PC3"), dim.to.compress=c("PC1", "PC2"), num.intervals= 20) {
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
