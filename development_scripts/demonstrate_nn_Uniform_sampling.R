library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)
library(FNN)
library(coda)

seed.number = 43
set.seed(seed.number)
dimensions <- c("PC1", "PC2", "PC3", "PC4", "PC5")
n.sample.points = 1000

# load data
data.dir <- "/home/mathis/Desktop/semesterarbeit10/data"

number.of.dimensions <- 15
env.data.raster.all.dim <- geodata::worldclim_global(var='bio', res=10, path=data.dir)  %>%
  terra::crop(terra::ext(-12, 25, 36, 60))

selected.layers <- sample(terra::nlyr(env.data.raster.all.dim), number.of.dimensions)
env.data.raster <- env.data.raster.all.dim[[selected.layers]]


# env.data.raster <- USE.MCMC::Worldclim_tmp %>%
#   terra::rast( type="xyz")

env.data.sf <- env.data.raster %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))

# Generate the environmental space using PCA
rpc <- rastPCA(env.data.raster,  stand = TRUE)
env.data.raster.with.pc <- c(env.data.raster, rpc$PCs)

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
virtual.presence.points.pc <- terra::extract(env.data.raster.with.pc, virtual.presence.points, bind = TRUE) %>%
  sf::st_as_sf()

# max resolution to expect uniform results. For resolutions higher we expect oversampling of dense regions
# lower resolutions should not be an issue, as oversampling of the border areas should be counteracted by removing points that had a
# distance larger than half a grid cell to their originally sampled point

print(paste("The resolution of the grid should be lower than ",maxResNn(env.data.raster = rpc,
               dimensions = dimensions,
               n.neighbors = 10,
               low.end.of.inclueded.points = 100, high.end.of.included.points = 4, PCA = TRUE)))

sampled.points.uniform.p <- paSampling(env.rast = env.data.raster, pres = virtual.presence.points, grid.res = 10)
sampled.points.uniform.p.location <- sf::st_drop_geometry(sampled.points.uniform.p) %>%
  dplyr::select(c("x", "y")) %>%
  sf::st_as_sf(coords = c("x", "y"))
sampled.points.uniform.paper <- terra::extract(env.data.raster.with.pc, sampled.points.uniform.p.location)

sampled.points.uniform.nn <- paSamplingNn(env.rast = env.data.raster,
                                          pres = virtual.presence.points,
                                          grid.res = 15, n.tr = 2,
                                          n.samples = n.sample.points)

sampled.points.mcmc <- paSamplingMcmc(env.data.raster = env.data.raster,
                                      pres = virtual.presence.points, precomputed.pca = rpc, environmental.cutof.percentile = 0.001,
                                      num.chains = 1,
                                      num.cores = 1,
                                      chain.length = 10000,
                                      dimensions = dimensions,
                                      n.samples = n.sample.points)




if(TRUE) {
  # Create a layout with 3 rows: 2 for plots, 1 for legend
  plotting.dimensions <- c("PC1", "PC2" , "PC3", "PC4", "PC5")
  layout(matrix(c(1:(length(plotting.dimensions) + 1)), nrow = length(plotting.dimensions) + 1, byrow = TRUE),
         heights = c(rep(2, length(plotting.dimensions)) ,1))  # Make legend row shorter

  par(mar = c(4, 4, 5, 2))  # bottom, left, top, right

  plots <- invisible(lapply(plotting.dimensions, function(col) {
    x_range <- range(
      c(density(env.with.pc.fs[[col]])$x,
        density(sampled.points.uniform.nn[[col]])$x,
        density(sampled.points.mcmc[[col]])$x)
    )
    y_range <- range(
      c(density(env.with.pc.fs[[col]])$y,
        density(sampled.points.uniform.nn[[col]])$y,
        density(sampled.points.mcmc[[col]])$y)
    )

    plot(density(env.with.pc.fs[[col]]),
         col = "green",
         main = paste("Density Comparison for", col),
         xlim = x_range,
         ylim = y_range,
         lwd = 2)

    lines(density(virtual.presence.points.pc[[col]]), col = "black", lwd = 2)
    lines(density(sampled.points.uniform.paper[[col]]), col = "orange", lwd = 2)
    lines(density(sampled.points.uniform.nn[[col]]), col = "blue", lwd = 2)
    lines(density(sampled.points.mcmc[[col]]), col = "red", lwd = 2)
  }))

  # Create an empty plot for the legend
  par(mar = c(0, 0, 0, 0))  # Remove margins
  plot.new()
  legend("center",
         legend = c("Environment", "Virtual Presence", "Sampled Points uniform",
                    "Sampled Points uniform nn", "Sampled Points mcmc"),
         col = c("green", "black", "orange", "blue", "red"),
         lwd = 2,
         horiz = FALSE,  # Vertical legend (rows)
         bty = "n",      # No box around legend
         cex = 1.5)      # incease text size
}


 # x <- analyseAreaInHigherDim(sf::st_drop_geometry(env.with.pc.fs))

library(ggplot2)
library(gridExtra)
library(cowplot)  # For get_legend() function

plotting.dimensions <- c("PC1", "PC2", "PC3", "PC4", "PC5")
line.with <- 0.8

# Create individual plots for each PC dimension
individual_plots <- lapply(plotting.dimensions, function(col) {
  ggplot() +
    geom_density(aes(x = env.with.pc.fs[[col]]), color = "green", size = line.with) +
    geom_density(aes(x = virtual.presence.points.pc[[col]]), color = "black", size = line.with) +
    geom_density(aes(x = sampled.points.uniform.paper[[col]]), color = "orange", size = line.with) +
    geom_density(aes(x = sampled.points.uniform.nn[[col]]), color = "blue", size = line.with) +
    geom_density(aes(x = sampled.points.mcmc[[col]]), color = "red", size = line.with) +
    labs(
      title = paste("Density Comparison for", col),
      x = "Value",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.position = "none"
    )
})

# Create shared legend
legend <- get_legend(
  plot <- ggplot() +
    geom_density(aes(x = c(1,2,3), color = "Environment"), size = line.with) +
    geom_density(aes(x = c(1,2,3), color = "Virtual Presence"), size = line.with) +
    geom_density(aes(x = c(1,2,3), color = "Sampled Points Uniform"), size = line.with) +
    geom_density(aes(x = c(1,2,3), color = "Sampled Points Uniform NN"), size = line.with) +
    geom_density(aes(x = c(1,2,3), color = "Sampled Points MCMC"), size = line.with) +
    scale_color_manual(
      name = "",
      values = c("Environment" = "green",
                 "Virtual Presence" = "black",
                 "Sampled Points Uniform" = "orange",
                 "Sampled Points Uniform NN" = "blue",
                 "Sampled Points MCMC" = "red")
    ) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12),
          legend.title = element_blank())
)

# Arrange plots in grid with legend at bottom
grid.arrange(
  arrangeGrob(grobs = individual_plots, ncol = 2),
  legend,
  heights = c(0.85, 0.15)
)


