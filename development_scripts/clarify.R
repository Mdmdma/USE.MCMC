library(ggplot2)
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(parallel)
library(FNN)
library(coda)

densityFunction <- function(point){
  if (point[[1]] < 5 && point[[1]] > 0){
    if (point[[2]] < 1 && point[[2]] > 0){
      return(1)
    }
  }
  return(0)
}

dimensions <- c("PC1", "PC2")
point = data.frame(PC1 = 0.5, PC2 = 0.5)

covariance.scaling <-0.1
covariance.matrix <- covariance.scaling * diag(1, nrow = length(dimensions))
proposalFunction <- addHighDimGaussian(cov.mat = covariance.scaling * covariance.matrix,
                                       dim = length(dimensions))
proposed.point <- proposalFunction(point = point, covariance.adjuster = 1, dim = dimensions)

sampled.points <- mcmcSampling(dataset = point,
                               dimensions = dimensions,
                               n.sample.points = 10000,
                               proposalFunction = proposalFunction,
                               densityFunction = densityFunction,
                               burnIn = TRUE,
                               verbose = TRUE)
hist(sampled.points[[1]])
hist(sampled.points[[2]])
