library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)
library(FNN)
library(coda)

seed.number = 51
dimensions <- c("PC1", "PC2")
set.seed(seed.number)

# load data
data.dir <- "/home/mathis/Desktop/semesterarbeit10/data"

number.of.dimensions <- 5
env.data.raster.all.dim <- geodata::worldclim_global(var='bio', res=10, path=data.dir)  %>%
  terra::crop(terra::ext(-12, 25, 36, 60))

selected.layers <- sample(terra::nlyr(env.data.raster.all.dim), number.of.dimensions)
env.data.raster <- env.data.raster.all.dim[[selected.layers]]


# env.data.raster <- USE.MCMC::Worldclim_tmp %>%
#   terra::rast( type="xyz")

env.data.sf <- env.data.raster %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))

# fixing the


# Generate the environmental space using PCA
rpc <- rastPCA(env.data.raster,  stand = TRUE)
raster.with.pc <- c(env.data.raster, rpc$PCs)

# Attaching the data in the PCA coordinates
env.with.pc.fs <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(env.data.sf)

# subsample env space to speed up the process
env.with.pc.fs.subsampled <- env.with.pc.fs[runif(min(nrow(env.with.pc.fs), 2000) , 1, nrow(env.with.pc.fs)),]
#specify the dimension that should be included in the following analysys

# Generate virtual species
virtual.presence.data <- getVirtualSpeciesPresencePoints(env.data = env.data.raster, n.samples = 300)
virtual.presence.points <- virtual.presence.data$sample.points
virtual.presence.points.pc <- terra::extract(raster.with.pc, virtual.presence.points, bind = TRUE) %>%
  sf::st_as_sf()

# max resolution to expect uniform results. For resolutions higher we expect oversampling of dense regions
# lower resolutions should not be an issue, as oversampling of the border areas should be counteracted by removing points that had a
# distance larger than half a grid cell to their originally sampled point

print(paste("The resolution of the grid should be lower than ",maxResNn(env.data.raster = rpc,
               dimensions = dimensions,
               n.neighbors = 10,
               low.end.of.inclueded.points = 100, high.end.of.included.points = 4, PCA = TRUE)))

sampled.points.uniform <- paSamplingNn(env.rast = env.data.raster, pres = virtual.presence.points, PCA = rpc)
sampled.points.mcmc <- paSamplingMcmc(env.data.raster = env.data.raster,
                                      pres = virtual.presence.points, precomputed.pca = rpc, environmental.cutof.percentile = 0.001)




if(TRUE) {
  # ggplot(sampled.points, aes(x = PC0, y = PC2)) +
  # geom_point() +
  # scale_color_viridis_c() +  # Optional: prettier color scale
  # theme_minimal()

  par(mfrow = c(1,2))
  plots <- invisible(lapply(dimensions, function(col) {
    # Create an empty plot with appropriate limits

    x_range <- range(
      c(density(env.with.pc.fs[[col]])$x,
        density(virtual.presence.points.pc[[col]])$x,
        density(sampled.points.uniform[[col]])$x,
        density(sampled.points.mcmc[[col]])$x)
    )
    y_range <- range(
      c(density(env.with.pc.fs[[col]])$y,
        density(virtual.presence.points.pc[[col]])$y,
        density(sampled.points.uniform[[col]])$y,
        density(sampled.points.mcmc[[col]])$y)
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

    # Add the denisty for the uniformly sampled points
    lines(density(sampled.points.uniform[[col]]),
          col = "orange",
          lwd = 2)

    # Add the denisty for the mcmc sampled points
    lines(density(sampled.points.mcmc[[col]]),
          col = "red",
          lwd = 2)

    # Add a legend
    legend("topleft",
           legend = c("Environment", "Virtual Presence", "Sampled Points uniform", "Sampled Points mcmc"),
           col = c("green", "black", "orange", "red"),
           lwd = 2)
  }))
}
