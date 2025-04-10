# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)


# Parameters:

#specify the dimension that should be included in the following analysys
dimensions <- c("PC1", "PC2") #, "PC3", "PC4","PC5"

#Nedded for plotting
par(mfrow = c(1, 1))

# load data
env.data.raster <- USE.MCMC::Worldclim_tmp %>%
  terra::rast( type="xyz")

# convert to SF dataframe
env.data.sf <- env.data.raster %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))

# fixing the
set.seed(42)


# Generate the environmental space using PCA
rpc <- rastPCA(env.data.raster,  stand = TRUE)

# Attaching the data in the PCA coordinates
env.with.pc.fs <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(env.data.sf)

# subsample env space to speed up the process
env.with.pc.fs <- env.with.pc.fs[runif(nrow(env.with.pc.fs)/10, 1, nrow(env.with.pc.fs)),]

# clean data
env.data.cleaned <- sf::st_drop_geometry(env.with.pc.fs[dimensions])


# environment model
environmental.data.model <- mclust::densityMclust(env.data.cleaned, plot = TRUE)
summary(environmental.data.model)
environmental.densities <- mclust::predict.densityMclust(environmental.data.model, env.data.cleaned)
environment.threshold <- stats::quantile(environmental.densities, 0.01)

# sample species model
virtual.presence.data <- getVirtualSpeciesPresencePoints(env.data = env.data.raster, n.samples = 300)
virtual.presence.points <- virtual.presence.data$sample.points
virtual.presence.points.pc <- terra::extract(rpc$PCs, virtual.presence.points, bind = TRUE) %>%
  sf::st_as_sf()

species.model = mclust::densityMclust(sf::st_drop_geometry(virtual.presence.points.pc[dimensions]), plot = TRUE)
summary(species.model)

#density Function
densityFunction <- mclustDensityFunction(env.model = environmental.data.model, presence.model = species.model,
                                         dim = dimensions, threshold = environment.threshold)

# set sampling parameters
covariance.proposal.function <-0.5
proposalFunction <- addHighDimGaussian(cov_mat =covariance.proposal.function * diag(length(dimensions)), dim = length(dimensions))


# sample points
sampled.points <- mcmcSampling(dataset = env.with.pc.fs[dimensions], dimensions = dimensions, n.sample.points = 1000,
                               proposalFunction = proposalFunction, densityFunction = densityFunction)


real.sampled.points <- sampled.points
real.sampled.points$PC3 <- 0
sampled.points$PC3 <- 0

#plot
par(mfrow = c(2, 2))
plot_points_with_lines(sampled.points, c("PC1", "PC2", "PC3"), title = paste("sampled points with coveariance", covariance.proposal.function),
                       limits = list(c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)), c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2))))
plot_points_with_lines(real.sampled.points, c("PC1", "PC2", "PC3"), title = paste("Covariance is diagonal ", covariance.proposal.function),
                       limits = list(c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)), c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2))))
plot(env.with.pc.fs$PC1, env.with.pc.fs$PC2, main = "Environment")
plot(virtual.presence.points.pc$PC1, virtual.presence.points.pc$PC2, xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
     ylim = c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)), main = " Virtual prescence points" )

par(mfrow = c(3, length(dimensions)))
invisible(lapply(dimensions, function(col) hist(env.with.pc.fs[[col]], main=paste("Histogram of environment", col))))
invisible(lapply(dimensions, function(col) hist(virtual.presence.points.pc[[col]], main=paste("Histogram of virtual presence", col))))
invisible(lapply(dimensions, function(col) hist(real.sampled.points[[col]], main=paste("Histogram of sampled points", col))))
par(mfrow = c(3, length(dimensions)))
invisible(lapply(dimensions, function(col) plot(density(env.with.pc.fs[[col]]), main=paste("Density of environment", col))))
invisible(lapply(dimensions, function(col) plot(density(virtual.presence.points.pc[[col]]), main=paste("Density of virtual presence", col))))
invisible(lapply(dimensions, function(col) plot(density(real.sampled.points[[col]]), main=paste("Density of sampled points", col))))
par(mfrow = c(1, 1))

#plotDensity2dpro(dataset =  real.sampled.points, species = virtual.presence.points.pc, xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)), ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
#                 densityFunction = densityFunction, resolution = 100)

#plotInGeographicalSpace1(presence.distribution.raster =  virtual.presence.data$original.distribution.raster,
#                         presence.points = virtual.presence.points.pc, absence.points = real.sampled.points )
