# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)


#Nedded for plotting
par(mfrow = c(1, 1))
plot <- TRUE

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

#specify the dimension that should be included in the following analysys
dimensions <- c("PC1", "PC2", "PC3") #, "PC3", "PC4","PC5"

# clean data
env.data.cleaned <- sf::st_drop_geometry(env.with.pc.fs[dimensions])


# environment model
environmental.data.model <- mclust::densityMclust(env.data.cleaned, plot = plot)
summary(environmental.data.model)
environmental.densities <- mclust::predict.densityMclust(environmental.data.model, env.data.cleaned)
environment.threshold <- stats::quantile(environmental.densities, 0.01)

# sample species model
virtual.presence.data <- getVirtualSpeciesPresencePoints(env.data = env.data.raster, n.samples = 300)
virtual.presence.points <- virtual.presence.data$sample.points
virtual.presence.points.pc <- terra::extract(rpc$PCs, virtual.presence.points, bind = TRUE) %>%
  sf::st_as_sf()

species.model = mclust::densityMclust(sf::st_drop_geometry(virtual.presence.points.pc[dimensions]),
                                      plot = plot)
summary(species.model)

#density Function
densityFunction <- mclustDensityFunction(env.model = environmental.data.model,
                                         presence.model = species.model,
                                         dim = dimensions,
                                         threshold = environment.threshold)

# # set sampling parameters
covariance.proposal.function <-0.075
proposalFunction <- addHighDimGaussian(cov.mat =covariance.proposal.function * diag(length(dimensions)),
                                       dim = length(dimensions))
#

# sample points
sampled.points <- mcmcSampling(dataset = env.with.pc.fs,
                               dimensions = dimensions,
                               n.sample.points = 10000,
                               proposalFunction = proposalFunction,
                               densityFunction = densityFunction)

# setup environment to compute in parallel
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
max.distance <- 2
clusterExport(cl, c("mapBackOnRealPoints", "env.with.pc.fs", "dimensions", "max.distance"))

# map back onto real points
real.sampled.points.list <- parallel::parApply(cl, sampled.points, 1, function(point) mapBackOnRealPoints(env.with.pc.fs, point, dimensions, threshold = max.distance))
stopCluster(cl)
cat(sum(is.na(real.sampled.points.list)) ,"points have no real counterpart in the environment space, given a maximal distance of ", max.distance, "!")
real.sampled.points.list.clean <- real.sampled.points.list[!is.na(real.sampled.points.list)]
real.sampled.points <- do.call(rbind, real.sampled.points.list.clean)
#plot
if (plot){
  par(mfrow = c(2, 2))
  plotPointsWithLines(sampled.points, c("PC1", "PC2", "PC3"),
                      title = paste("sampled points with coveariance", covariance.proposal.function),
                      limits = list(c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)), c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2))))
  plotPointsWithLines(real.sampled.points, c("PC1", "PC2", "PC3"),
                      title = paste("Covariance is diagonal ", covariance.proposal.function),
                      limits = list(c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)), c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2))))
  plot(env.with.pc.fs$PC1, env.with.pc.fs$PC2, main = "Environment")
  plot(virtual.presence.points.pc$PC1, virtual.presence.points.pc$PC2,
       xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
       ylim = c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
       main = " Virtual prescence points" )

  par(mfrow = c(3, length(dimensions)))
  invisible(lapply(dimensions, function(col) hist(env.with.pc.fs[[col]],
                                                  main=paste("Histogram of environment", col))))
  invisible(lapply(dimensions, function(col) hist(virtual.presence.points.pc[[col]],
                                                  main=paste("Histogram of virtual presence", col))))
  invisible(lapply(dimensions, function(col) hist(real.sampled.points[[col]],
                                                  main=paste("Histogram of sampled points", col))))
  par(mfrow = c(3, length(dimensions)))
  invisible(lapply(dimensions, function(col) plot(density(env.with.pc.fs[[col]]),
                                                  main=paste("Density of environment", col))))
  invisible(lapply(dimensions, function(col) plot(density(virtual.presence.points.pc[[col]]),
                                                  main=paste("Density of virtual presence", col))))
  invisible(lapply(dimensions, function(col) plot(density(real.sampled.points[[col]]), main=paste("Density of sampled points", col))))
  par(mfrow = c(1, 1))

  plotDensity2dpro(dataset =  real.sampled.points,
                   species = virtual.presence.points.pc,
                   xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                   ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                   densityFunction = densityFunction,
                   resolution = 100)

  plotInGeographicalSpace(presence.distribution.raster =  virtual.presence.data$original.distribution.raster,
                          presence.points = virtual.presence.points.pc,
                          absence.points = real.sampled.points )
}
