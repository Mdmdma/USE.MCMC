acceptNextPoint <- function(current.point, proposed.point, densityFunction){
  acceptance.ratio <- densityFunction(proposed.point) / densityFunction(current.point)
  if (stats::runif(1) < acceptance.ratio) return(TRUE)
  return(FALSE)
}
