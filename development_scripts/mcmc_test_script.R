# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)
library(FNN)


#Needed for plotting
par(mfrow = c(1, 1))
plot <- TRUE

# load data
#env.data.raster <- USE.MCMC::Worldclim_tmp %>%
#  terra::rast( type="xyz")
datadir <- "/home/mathis/Desktop/semesterarbeit10/data"
env.data.raster <- geodata::worldclim_global(var='bio', res=2.5, path=datadir)  %>%
  terra::crop(terra::ext(-12, 25, 36, 60))
  #terra::rast( type="xyz")

env.data.raster <- geodata::worldclim_country(country = "ch", var = "bio", path=datadir, res=0.5)

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
dimensions <- c("PC1", "PC2", "PC3", "PC4","PC5") #, "PC3", "PC4","PC5"

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

mapped.sampled.point.locations <- FNN::get.knnx(env.data.cleaned[dimensions], sampled.points[dimensions],k = 1)
mapped.sampled.points <- env.with.pc.fs[mapped.sampled.point.locations$nn.index,]
mapped.sampled.points$density <- sampled.points$density
mapped.sampled.points$distance <- mapped.sampled.point.locations$nn.dist

distance.threshold <- stats::quantile(mapped.sampled.points$distance, 0.95)
filtered.mapped.sampled.points <- mapped.sampled.points[mapped.sampled.points$distance < distance.threshold, ]


n.samples <- 3000
sample.indexes <- indices <- seq(1, nrow(filtered.mapped.sampled.points), length.out = n.samples)
real.sampled.points <- filtered.mapped.sampled.points[indices, ]
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

  par(mfrow = c(length(dimensions),1))
  invisible(lapply(dimensions, function(col) {
    # Create an empty plot with appropriate limits
    x_range <- range(
      c(density(env.with.pc.fs[[col]])$x,
        density(virtual.presence.points.pc[[col]])$x,
        density(real.sampled.points[[col]])$x)
    )
    y_range <- range(
      c(density(env.with.pc.fs[[col]])$y,
        density(virtual.presence.points.pc[[col]])$y,
        density(real.sampled.points[[col]])$y)
    )

    # Plot the first density
    plot(density(env.with.pc.fs[[col]]),
         col = "green",
         main = paste("Density Comparison for", col),
         xlim = x_range,
         ylim = y_range,
         lwd = 2)

    # Add the second density to the same plot
    lines(density(virtual.presence.points.pc[[col]]),
          col = "black",
          lwd = 2)

    # Add the third density to the same plot
    lines(density(real.sampled.points[[col]]),
          col = "red",
          lwd = 2)

    # Add a legend
    legend("topleft",
           legend = c("Environment", "Virtual Presence", "Sampled Points"),
           col = c("green", "black", "red"),
           lwd = 2)
  }))
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

