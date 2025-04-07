# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)


#Nedded for plotting
par(mfrow = c(1, 1))

# load data
envDataRaster <- USE.MCMC::Worldclim_tmp %>%
  terra::rast( type="xyz")

envDataSf <- envDataRaster %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))

# Create virtual species
set.seed(153)


# Generate the environmental space using PCA
rpc <- rastPCA(envDataRaster,  stand = TRUE)


envWithPcSf <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(envDataSf)

# subsample env space to speed up the process
envWithPcSf <- envWithPcSf[runif(nrow(envWithPcSf)/10, 1, nrow(envWithPcSf)),]

dimensions <- c("PC1", "PC2") #, "PC3", "PC4","PC5"
# cleaned data
environmentalData <- sf::st_drop_geometry(envWithPcSf[dimensions])


# environment model
environmental.data.model <- mclust::densityMclust(environmentalData, plot = TRUE)
summary(environmental.data.model)
environmental.densities <- mclust::predict.densityMclust(environmental.data.model, environmentalData)
environment.threshold <- stats::quantile(environmental.densities, 0.01)

# sample species model
virtual.precence.points <- getVirtualSpeciesPresencePoints(environemtalData = envDataRaster, n.samples = 1000)
virtual.precence.points.pc <- terra::extract(rpc$PCs, virtual.precence.points, bind = TRUE) %>%
  sf::st_as_sf()
virtual.precence.points.pc <- sf::st_drop_geometry(virtual.precence.points.pc[dimensions])

species.model = mclust::densityMclust(virtual.precence.points.pc, plot = TRUE)
summary

#density Function
densityFunction <- mclustDensityFunction(environmentalModel = environmental.data.model, presenceModel = species.model,
                                         dim = dimensions, threshold = environment.threshold)

# set sampling parameters
covariance.proposal.function <-0.3
proposalFunction <- addHighDimGaussian(cov_mat =covariance.proposal.function * diag(length(dimensions)), dim = length(dimensions))


# sample points
sampled.points <- mcmcSampling(dataset = envWithPcSf, dimensions = dimensions, n.sample.points = 100,
                               proposalFunction = proposalFunction, densityFunction = densityFunction)


#plot

par(mfrow = c(3, 1))
plot_points_with_lines(sampled.points, c("PC1", "PC2", "PC3"), title = paste("Covariance is diagonal ", covariance.proposal.function),
                       limits = list(c(min(envWithPcSf$PC1), max(envWithPcSf$PC1)), c(min(envWithPcSf$PC2), max(envWithPcSf$PC2))))
plot(envWithPcSf$PC1, envWithPcSf$PC2, main = paste("Covariance is diagonal ", covariance.proposal.function))
plot(virtual.precence.points.pc$PC1, virtual.precence.points.pc$PC2, xlim = c(min(envWithPcSf$PC1), max(envWithPcSf$PC1)),
     ylim = c(min(envWithPcSf$PC2), max(envWithPcSf$PC2)) )
par(mfrow = c(2, length(dimensions)))
invisible(lapply(dimensions, function(col) hist(envWithPcSf[[col]], main=paste("Histogram of environment", col))))
invisible(lapply(dimensions, function(col) hist(sampled.points[[col]], main=paste("Histogram of sampled points", col))))


