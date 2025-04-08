# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)


#Nedded for plotting
par(mfrow = c(1, 1))

# load data
envDataRaster <- USE.MCMC::Worldclim_tmp %>%
  terra::rast( type="xyz")

# convert to SF dataframe
envDataSf <- envDataRaster %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))

# Create virtual species
set.seed(42)


# Generate the environmental space using PCA
rpc <- rastPCA(envDataRaster,  stand = TRUE)

# Attaching the data in the PCA coordinates
envWithPcSf <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(envDataSf)

# subsample env space to speed up the process
envWithPcSf <- envWithPcSf[runif(nrow(envWithPcSf)/10, 1, nrow(envWithPcSf)),]

#specify the dimension that should be included in the following analysys
dimensions <- c("PC1", "PC2") #, "PC3", "PC4","PC5"

# clean data
environmentalDataCleaned <- sf::st_drop_geometry(envWithPcSf[dimensions])


# environment model
environmental.data.model <- mclust::densityMclust(environmentalDataCleaned, plot = TRUE)
summary(environmental.data.model)
environmental.densities <- mclust::predict.densityMclust(environmental.data.model, environmentalDataCleaned)
environment.threshold <- stats::quantile(environmental.densities, 0.05)

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
covariance.proposal.function <-0.1
proposalFunction <- addHighDimGaussian(cov_mat =covariance.proposal.function * diag(length(dimensions)), dim = length(dimensions))


# sample points
sampled.points <- mcmcSampling(dataset = envWithPcSf, dimensions = dimensions, n.sample.points = 1000,
                               proposalFunction = proposalFunction, densityFunction = densityFunction)

# setup environment to compute in parallel
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
max.distance <- 2
clusterExport(cl, c("mapBackOnRealPoints", "envWithPcSf", "dimensions", "max.distance"))

# map back onto real points
real.sampled.points.list <- parallel::parApply(cl, sampled.points, 1, function(point) mapBackOnRealPoints(envWithPcSf, point, dimensions, threshold = max.distance))
cat(sum(is.na(real.sampled.points.list)) ,"points have no real counterpart in the environment space, given a maximal distance of ", max.distance, "!")
real.sampled.points.list.clean <- real.sampled.points.list[!is.na(real.sampled.points.list)]
real.sampled.points <- do.call(rbind, real.sampled.points.list.clean)
#plot
par(mfrow = c(3, 1))
plot_points_with_lines(real.sampled.points, c("PC1", "PC2", "PC3"), title = paste("Covariance is diagonal ", covariance.proposal.function),
                       limits = list(c(min(envWithPcSf$PC1), max(envWithPcSf$PC1)), c(min(envWithPcSf$PC2), max(envWithPcSf$PC2))))
plot(envWithPcSf$PC1, envWithPcSf$PC2, main = paste("Covariance is diagonal ", covariance.proposal.function))
plot(virtual.precence.points.pc$PC1, virtual.precence.points.pc$PC2, xlim = c(min(envWithPcSf$PC1), max(envWithPcSf$PC1)),
     ylim = c(min(envWithPcSf$PC2), max(envWithPcSf$PC2)) )
par(mfrow = c(3, length(dimensions)))
invisible(lapply(dimensions, function(col) hist(envWithPcSf[[col]], main=paste("Histogram of environment", col))))
invisible(lapply(dimensions, function(col) hist(virtual.precence.points.pc[[col]], main=paste("Histogram of virtual presence", col))))
invisible(lapply(dimensions, function(col) hist(real.sampled.points[[col]], main=paste("Histogram of sampled points", col))))
par(mfrow = c(3, length(dimensions)))
invisible(lapply(dimensions, function(col) plot(density(envWithPcSf[[col]]), main=paste("Histogram of environment", col))))
invisible(lapply(dimensions, function(col) plot(density(virtual.precence.points.pc[[col]]), main=paste("Histogram of virtual presence", col))))
invisible(lapply(dimensions, function(col) plot(density(real.sampled.points[[col]]), main=paste("Histogram of sampled points", col))))


