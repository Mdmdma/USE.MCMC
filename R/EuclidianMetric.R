euclidianMetric <- function(pointA=NULL, pointB =NULL, dim = ""){
  pointA <- sf::st_drop_geometry(pointA)
  pointB <- sf::st_drop_geometry(pointB)
  m <- Map(`-`, pointA[dim], pointB[dim]) %>% unlist() %>% sum() %>% abs()
}
