mapBackOnRealPoints <- function(dataset = NULL, point = NULL, dim ="", threshold = 0.5){
  distances <- apply(dataset, 1, function(row) euclidianMetric(row, point, dim))
  min.dist.index <- which.min(distances)
  min.dist <- distances[min.dist.index]
  best.point <- dataset[min.dist.index,]
  best.point$distanceFromSampledPoint <- min.dist
  if (min.dist < threshold) return(best.point)
  return(NA)
}
