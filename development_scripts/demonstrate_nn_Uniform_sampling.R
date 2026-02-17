library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)
library(FNN)
library(coda)
library(profvis)

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
env.with.pc.sf <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(env.data.sf)

# subsample env space to speed up the process
env.with.pc.sf.subsampled <- env.with.pc.sf[runif(min(nrow(env.with.pc.sf), 2000) , 1, nrow(env.with.pc.sf)),]
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

to.profile <- function(nn.based.exclustion = TRUE){
  if (nn.based.exclustion){
    optimalDistanceThresholdNn(env.data = env.with.pc.sf)
  }
  paSamplingNn(env.rast = env.data.raster,
               pres = virtual.presence.points,
               grid.res = 15, n.tr = 2,
               n.samples = n.sample.points,
               #precomputed.pca = rpc,
               nn.based.presence.exclusion = nn.based.exclustion)
}


profvis({
  for (i in 1:5) {
    to.profile(FALSE)
    }
  }
)

profvis({
  for (i in 1:1) {
    x <- paSampling(env.rast = env.data.raster,
                 pres = virtual.presence.points,
                 grid.res = 15, n.tr = 2)
  }
}
)


sampled.points.uniform.nn <- paSamplingNn(env.rast = env.data.raster,
                                          pres = virtual.presence.points,
                                          grid.res = 15, n.tr = 2,
                                          n.samples = n.sample.points)




sampled.points.mcmc <- paSamplingMcmc(env.data.raster = env.data.raster,
                                      pres = virtual.presence.points, precomputed.pca = rpc, environmental.cutof.percentile = 0.001,
                                      num.chains = 1,
                                      num.cores = 1,
                                      covariance.correction = 35,
                                      chain.length = 10000,
                                      dimensions = dimensions,
                                      n.samples = n.sample.points)


env.with.pc.sf.subsampled <- env.with.pc.sf[
  stats::runif(min(nrow(env.with.pc.sf), 2000) , 1, nrow(env.with.pc.sf)),]
env.data.cleaned.subsampled <- sf::st_drop_geometry(
  env.with.pc.sf.subsampled[dimensions])
environmental.data.model <- mclust::densityMclust(env.data.cleaned.subsampled)
title("Gaussian mixture fitted to the environment")

selected.dimensions <- c("PC1", "PC2")
env.with.pc <- env.with.pc.sf %>%
  sf::st_drop_geometry()
test.point <- env.with.pc[1, ]

profvis({
  for (i in 1:1000) {
    x <- get.nearest.ten.neighbors <- FNN::get.knnx(env.with.pc[selected.dimensions],
                                               test.point[selected.dimensions],
                                               k =100)
  }
}
)



num.neighbors <- 10
local.nn.based.local.environmental.cutoff <- function(env.with.pc, num.neighbors){
  at.point.estimator <- function(point){
    nearest.neighbors.locations <- FNN::get.knnx(env.with.pc[selected.dimensions],
                                           point[selected.dimensions],
                                           k=2)
    mean.dist.to.neighbors <- mean(nearest.neighbors.locations$nn.dist)
    print(paste("mean dist to neighbors ", mean.dist.to.neighbors))
    nearest.neighbors <- env.with.pc[nearest.neighbors.locations$nn.index,]
    stats.of.nn <- FNN::get.knnx(env.with.pc[selected.dimensions],
                                 nearest.neighbors[selected.dimensions],
                                 k=num.neighbors)
    mean.dist.to.neighbors.in.neighborhood <- mean(stats.of.nn$nn.dist[1:2,])
    print(paste("mean dist in neighbrohood ", mean.dist.to.neighbors.in.neighborhood))
    if (mean.dist.to.neighbors <= mean.dist.to.neighbors.in.neighborhood * 1.2){
      output <- 1
    } else {
      output <- 0
    }
  }
}

density.function <- local.nn.based.local.environmental.cutoff(env.with.pc = env.with.pc, num.neighbors = 100)
x <- density.function(test.point)

gridsize <- 30
grid <- expand.grid(
  x = seq(min(env.with.pc$PC1), max(env.with.pc$PC1), length.out = gridsize),
  y = seq(min(env.with.pc$PC2), max(env.with.pc$PC2), length.out = gridsize)
)
names(grid) <- c("PC1", "PC2")

grid$density <- purrr::map_dbl(1:nrow(grid), ~ density.function(grid[.x, , drop = FALSE]))

p.density.naive <- ggplot() +
  geom_tile(data=grid, aes(x=PC1, y=PC2, fill=density)) +
  geom_point(data=env.with.pc, aes(x=PC1, y=PC2), size = 0.1, alpha=0.2) +
  theme_minimal() +
  labs(title = "Density function in the naive case")
print(p.density.naive)






