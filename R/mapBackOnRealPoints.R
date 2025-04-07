mapBackOnRealPoints <- function(dataset = NULL, point = NULL, dim ="", threshold = 0.1){
  distances <- apply(dataset, 1, function(row) euclidianMetric(row, point, dim))
  min.dist.index <- which.min(distances)
  min.dist <- distances[min.dist.index]
  current.point <- dataset[min.dist.index,]
  if (min.dist < threshold) return(dataset[min.dist.index,])
  return(NA)
}
