mcmcSampling <- function(dataset = NULL, dimensions= list(""), densityFunction = alwaysOne, proposalFunction = addHighDimGaussian , n.sample.points = 0){
  starting.index <- stats::runif(1,1,nrow(dataset))
  current.point <- dataset[starting.index,]
  dataset <- dataset[-starting.index,]
  sampled.points <- dataset[0, ]
  while (nrow(sampled.points) < n.sample.points) {
    print(paste("sampling point number", nrow(sampled.points)))
    proposed.point <- addHighDimGaussian(currentPoint = current.point, dim = dimensions)
    if (acceptNextPoint(current.point, proposed.point, densityFunction)){
      distances <- apply(dataset, 1, function(row) euclidianMetric(row, current.point, dimensions))
      min.dist.index <- which.min(distances)
      current.point <- dataset[min.dist.index,]
      dataset <- dataset[-min.dist.index,]
      sampled.points <- rbind(sampled.points, current.point)
      print(nrow(dataset))
    }
  }
  return(sampled.points)
}
