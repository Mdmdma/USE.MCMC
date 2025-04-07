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
set.seed(153)


# Generate the environmental space using PCA
rpc <- rastPCA(envData,  stand = TRUE)


envWithPc <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(sfEnvData)

# subsample env space to speed up the process
envWithPc <- envWithPc[runif(nrow(envWithPc)/10, 1, nrow(envWithPc)),]

dimensions <- c("PC1", "PC2") #, "PC3", "PC4","PC5"
# cleaned data
environmentalData <- sf::st_drop_geometry(envWithPc[dimensions])


# environment model
environmental.data.model <- mclust::densityMclust(environmentalData, plot = TRUE)
summary(environmental.data.model)
environmental.densities <- mclust::predict.densityMclust(environmental.data.model, environmentalData)
threshold <- stats::quantile(environmental.densities, 0.01)

# sample species model
virtual.precence.points <- getVirtualSpeciesPresencePoints(environemtalData = envData, n.samples = 1000)
virtual.precence.points.pc <- terra::extract(rpc$PCs, virtual.precence.points, bind = TRUE) %>%
  sf::st_as_sf()
virtual.precence.points.pc <- sf::st_drop_geometry(virtual.precence.points.pc[dimensions])

species.model = mclust::densityMclust(virtual.precence.points.pc, plot = TRUE)
summary

#density Function
densityFunction <- mclustDensityFunction(environmentalModel = environmental.data.model, presenceModel = species.model,
                                         dim = dimensions, threshold = threshold)

# set sampling parameters
covariance <-0.3
proposalFunction <- addHighDimGaussian(cov_mat =covariance * diag(length(dimensions)), dim = length(dimensions))


# sample points
sampled.points <- mcmcSampling(dataset = envWithPc, dimensions = dimensions, n.sample.points = 1000,
                               proposalFunction = proposalFunction, densityFunction = densityFunction)


#plot

par(mfrow = c(3, 1))
plot_points_with_lines(sampled.points, c("PC1", "PC2", "PC3"), title = paste("Covariance is diagonal ", covariance),
                       limits = list(c(min(envWithPc$PC1), max(envWithPc$PC1)), c(min(envWithPc$PC2), max(envWithPc$PC2))))
plot(envWithPc$PC1, envWithPc$PC2, main = paste("Covariance is diagonal ", covariance))
plot(virtual.precence.points.pc$PC1, virtual.precence.points.pc$PC2, xlim = c(min(envWithPc$PC1), max(envWithPc$PC1)),
     ylim = c(min(envWithPc$PC2), max(envWithPc$PC2)) )
par(mfrow = c(2, length(dimensions)))
invisible(lapply(dimensions, function(col) hist(envWithPc[[col]], main=paste("Histogram of environment", col))))
invisible(lapply(dimensions, function(col) hist(sampled.points[[col]], main=paste("Histogram of sampled points", col))))


