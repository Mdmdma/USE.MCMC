# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)

#Nedded for plotting
par(mfrow = c(1, 1))

# load data
envData <- USE.MCMC::Worldclim_tmp %>%
  terra::rast( type="xyz")

sfEnvData <- envData %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))

# Create virtual species
set.seed(123)


# Generate the environmental space using PCA
rpc <- rastPCA(envData,  stand = TRUE)


envWithPc <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(sfEnvData)

# subsample env space to speed up the process
envWithPc <- envWithPc[runif(nrow(envWithPc)/10, 1, nrow(envWithPc)),]

# sample points
dimensions <- c("PC1", "PC2", "PC3")
covariance <- 1
proposalFunction <- addHighDimGaussian(cov_mat =covariance * diag(length(dimensions)), dim = length(dimensions))
sampled.points <- mcmcSampling(dataset = envWithPc, dimensions = dimensions, n.sample.points = 100, proposalFunction = proposalFunction)

#plot

par(mfrow = c(2, 1))
plot_points_with_lines(sampled.points, c("PC1", "PC2"), title = paste("Covariance is diagonal ", covariance),
                       limits = list(c(min(envWithPc$PC1), max(envWithPc$PC1)), c(min(envWithPc$PC2), max(envWithPc$PC2))))
plot(envWithPc$PC1, envWithPc$PC2, main = paste("Covariance is diagonal ", covariance))
par(mfrow = c(2, 3))
hist(envWithPc$PC1)
hist(envWithPc$PC2)
hist(envWithPc$PC3)
hist(sampled.points$PC1)
hist(sampled.points$PC2)
hist(sampled.points$PC3)



