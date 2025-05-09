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
plot <- FALSE

datadir <- "/home/mathis/Desktop/semesterarbeit10/data"
# load data
env.data.raster <- USE.MCMC::Worldclim_tmp %>%
  terra::rast( type="xyz")

# env.data.raster <- geodata::worldclim_global(var='bio', res=10, path=datadir)  %>%
#   terra::crop(terra::ext(-12, 25, 36, 60))
  #terra::rast( type="xyz")
#
# env.data.raster <- geodata::worldclim_country(country = "ch", var = "bio", path=datadir, res=2.5)
# # convert to SF dataframe
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
env.with.pc.fs <- env.with.pc.fs[runif(min(nrow(env.with.pc.fs), 2000) , 1, nrow(env.with.pc.fs)),]

#specify the dimension that should be included in the following analysys
dimensions <- c("PC1", "PC2") #, "PC3", "PC4","PC5"

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
proposalFunction <- addHighDimGaussian(cov.mat =covariance.proposal.function * diag(rpc$pca$sdev[1:length(dimensions)]),
                                       dim = length(dimensions))
#

# sample points
sampled.points <- mcmcSampling(dataset = env.with.pc.fs,
                               dimensions = dimensions,
                               n.sample.points = 10000,
                               proposalFunction = proposalFunction,
                               densityFunction = densityFunction,
                               burnIn = TRUE)

mapped.sampled.point.locations <- FNN::get.knnx(env.data.cleaned[dimensions], sampled.points[dimensions],k = 1)
mapped.sampled.points <- env.with.pc.fs[mapped.sampled.point.locations$nn.index,]
mapped.sampled.points$density <- sampled.points$density
mapped.sampled.points$distance <- mapped.sampled.point.locations$nn.dist

distance.threshold <- stats::quantile(mapped.sampled.points$distance, 0.95)
filtered.mapped.sampled.points <- mapped.sampled.points[mapped.sampled.points$distance < distance.threshold, ]


n.samples <- 500
sample.indexes <- floor(seq(1, nrow(filtered.mapped.sampled.points), length.out = min(n.samples, nrow(filtered.mapped.sampled.points))))
real.sampled.points <- filtered.mapped.sampled.points[indices, ]
# TODO check where the coordinate system gets lost
st_crs(real.sampled.points) <- 4326
#plot
if (plot){
  l1 <- plotDensityLines(dataset = sampled.points,
                         lines = TRUE, cols = c("PC1", "PC2", "PC3"),
                         density = FALSE, densityFunction = densityFunction, resolution = 100, species = virtual.presence.points.pc,
                         xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                         ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                         title = NULL)

  l2 <- plotDensityLines(dataset = real.sampled.points,
                         lines = TRUE, cols = c("PC1", "PC2", "PC3"),
                         density = FALSE, densityFunction = densityFunction, resolution = 100, species = virtual.presence.points.pc,
                         xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                         ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                         title = NULL)

  env.scatterplot <- ggplot(env.with.pc.fs, aes(x = PC1, y = PC2)) +
    geom_point() +
    ggtitle("Environment") +
    theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5))

  species.scatterplot <- ggplot(virtual.presence.points.pc, aes(x = PC1, y = PC2)) +
    geom_point() +
    xlim(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)) +
    ylim(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)) +
    ggtitle("Virtual presence points") +
    theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5))

  cowplot::plot_grid(l1, l2, env.scatterplot, species.scatterplot, ncol = 2)
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

  plotDensityLines(dataset = real.sampled.points,
                   lines = FALSE, cols = c("PC1", "PC2", "PC3"),
                   density = TRUE, densityFunction = densityFunction, resolution = 100, species = virtual.presence.points.pc,
                   xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                   ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                   title = NULL)

  plotInGeographicalSpace(presence.distribution.raster =  virtual.presence.data$original.distribution.raster,
                          presence.points = virtual.presence.points.pc,
                          absence.points = real.sampled.points )
}

