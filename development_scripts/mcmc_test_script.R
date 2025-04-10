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
env.data.pc.sf <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(env.data.sf)

# subsample env space to speed up the process
env.data.pc.sf <- env.data.pc.sf[runif(nrow(env.data.pc.sf)/10, 1, nrow(env.data.pc.sf)),]

#specify the dimension that should be included in the following analysys
dimensions <- c("PC1", "PC2") #, "PC3", "PC4","PC5"

# clean data
env.data.cleaned <- sf::st_drop_geometry(env.data.pc.sf[dimensions])


# environment model
environmental.data.model <- mclust::densityMclust(env.data.cleaned, plot = TRUE)
summary(environmental.data.model)
environmental.densities <- mclust::predict.densityMclust(environmental.data.model, env.data.cleaned)
environment.threshold <- stats::quantile(environmental.densities, 0.05)

# sample species model
virtual.precence.points <- getVirtualSpeciesPresencePoints(env.data = env.data.raster, n.samples = 100)
virtual.precence.points.pc <- terra::extract(rpc$PCs, virtual.precence.points, bind = TRUE) %>%
  sf::st_as_sf()
virtual.precence.points.pc <- sf::st_drop_geometry(virtual.precence.points.pc[dimensions])

species.model = mclust::densityMclust(virtual.precence.points.pc, plot = TRUE)
summary

#density Function
densityFunction <- mclustDensityFunction(env.model = environmental.data.model, presence.model = species.model,
                                         dim = dimensions, threshold = environment.threshold)

# set sampling parameters
covariance.proposal.function <- 0.5
proposalFunction <- addHighDimGaussian(cov_mat =covariance.proposal.function * diag(length(dimensions)), dim = length(dimensions))


# sample points
sampled.points <- mcmcSampling(dataset = env.data.pc.sf, dimensions = dimensions, n.sample.points = 100,
                               proposalFunction = proposalFunction, densityFunction = densityFunction)

# setup environment to compute in parallel
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
max.distance <- 2
clusterExport(cl, c("mapBackOnRealPoints", "env.data.pc.sf", "dimensions", "max.distance"))

# map back onto real points
real.sampled.points.list <- parallel::parApply(cl, sampled.points, 1, function(point) mapBackOnRealPoints(env.data.pc.sf, point, dimensions, threshold = max.distance))
cat(sum(is.na(real.sampled.points.list)) ,"points have no real counterpart in the environment space, given a maximal distance of ", max.distance, "!")
real.sampled.points.list.clean <- real.sampled.points.list[!is.na(real.sampled.points.list)]
real.sampled.points <- do.call(rbind, real.sampled.points.list.clean)
#plot
par(mfrow = c(2, 2))
plot_points_with_lines(sampled.points, c("PC1", "PC2", "PC3"), title = paste("sampled points with coveariance", covariance.proposal.function),
                       limits = list(c(min(env.data.pc.sf$PC1), max(env.data.pc.sf$PC1)), c(min(env.data.pc.sf$PC2), max(env.data.pc.sf$PC2))))
plot_points_with_lines(real.sampled.points, c("PC1", "PC2", "PC3"), title = paste("Covariance is diagonal ", covariance.proposal.function),
                       limits = list(c(min(env.data.pc.sf$PC1), max(env.data.pc.sf$PC1)), c(min(env.data.pc.sf$PC2), max(env.data.pc.sf$PC2))))
plot(env.data.pc.sf$PC1, env.data.pc.sf$PC2, main = "Environment")
plot(virtual.precence.points.pc$PC1, virtual.precence.points.pc$PC2, xlim = c(min(env.data.pc.sf$PC1), max(env.data.pc.sf$PC1)),
     ylim = c(min(env.data.pc.sf$PC2), max(env.data.pc.sf$PC2)), main = " Virtual prescence points" )

par(mfrow = c(3, length(dimensions)))
invisible(lapply(dimensions, function(col) hist(env.data.pc.sf[[col]], main=paste("Histogram of environment", col))))
invisible(lapply(dimensions, function(col) hist(virtual.precence.points.pc[[col]], main=paste("Histogram of virtual presence", col))))
invisible(lapply(dimensions, function(col) hist(real.sampled.points[[col]], main=paste("Histogram of sampled points", col))))
par(mfrow = c(3, length(dimensions)))
invisible(lapply(dimensions, function(col) plot(density(env.data.pc.sf[[col]]), main=paste("Density of environment", col))))
invisible(lapply(dimensions, function(col) plot(density(virtual.precence.points.pc[[col]]), main=paste("Density of virtual presence", col))))
invisible(lapply(dimensions, function(col) plot(density(real.sampled.points[[col]]), main=paste("Density of sampled points", col))))
par(mfrow = c(1, 1))

plotDensity2dpro(dataset =  real.sampled.points, species = virtual.precence.points.pc, xlim = c(min(env.data.pc.sf$PC1), max(env.data.pc.sf$PC1)),
                 ylim =c(min(env.data.pc.sf$PC2), max(env.data.pc.sf$PC2)),
                 densityFunction = densityFunction, resolution = 100)

