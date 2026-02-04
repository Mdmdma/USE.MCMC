# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)
library(FNN)
library(coda)


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
env.with.pc.fs.subsampled <- env.with.pc.fs[runif(min(nrow(env.with.pc.fs), 2000) , 1, nrow(env.with.pc.fs)),]
#specify the dimension that should be included in the following analysys
dimensions <- c("PC1", "PC2") #, "PC3", "PC4","PC5"

# clean data
env.data.cleaned <- sf::st_drop_geometry(env.with.pc.fs[dimensions])
env.data.cleaned.subsampled <- sf::st_drop_geometry(env.with.pc.fs.subsampled[dimensions])


# environment model
environmental.data.model <- mclust::densityMclust( env.data.cleaned.subsampled, plot = plot)
summary(environmental.data.model)
environmental.densities <- mclust::predict.densityMclust(environmental.data.model,  env.data.cleaned.subsampled)
environment.threshold <- stats::quantile(environmental.densities, 0.01)

# sample species model
virtual.presence.data <- getVirtualSpeciesPresencePoints(env.data = env.data.raster, n.samples = 300)
virtual.presence.points <- virtual.presence.data$sample.points
virtual.presence.points.pc <- terra::extract(rpc$PCs, virtual.presence.points, bind = TRUE) %>%
  sf::st_as_sf()

species.model = mclust::densityMclust(sf::st_drop_geometry(virtual.presence.points.pc[dimensions]),
                                      plot = plot)
summary(species.model)

species.model = mclust::densityMclust(sf::st_drop_geometry(virtual.presence.points.pc[dimensions]))
summary(species.model)
species.densities <- species.model$density
species.cutoff.threshold <- stats::quantile(species.densities, 0.9)
densityFunction <- mclustDensityFunction(env.model = environmental.data.model,
                                         species.model = species.model,
                                         dim = dimensions,
                                         threshold = environment.threshold,
                                         species.cutoff.threshold = species.cutoff.threshold)

# # set sampling parameters
covariance.proposal.function <-0.075
proposalFunction <- addHighDimGaussian(cov.mat =covariance.proposal.function * diag(rpc$pca$sdev[1:length(dimensions)]),
                                       dim = length(dimensions))
#

# sample points
n_cores <- 20
num.chains <- 4
chain.list <- list()
chain.list <- mclapply(1:num.chains, function(i) {
  sampled.points <- mcmcSampling(dataset = env.with.pc.fs.subsampled,
                                dimensions = dimensions,
                                n.sample.points = 50000,
                                proposalFunction = proposalFunction,
                                densityFunction = densityFunction,
                                burnIn = 1001,
                                covariance.correction = 110)
  mapped.sampled.point.locations <- FNN::get.knnx(env.data.cleaned[dimensions], sampled.points[dimensions],k = 1)
  mapped.sampled.points <- env.with.pc.fs[mapped.sampled.point.locations$nn.index,]
  mapped.sampled.points$density <- sampled.points$density
  mapped.sampled.points$distance <- mapped.sampled.point.locations$nn.dist

  distance.threshold <- stats::quantile(mapped.sampled.points$distance, 0.95)
  filtered.mapped.sampled.points <- mapped.sampled.points[mapped.sampled.points$distance < distance.threshold, ]


  n.samples <- 500
  sample.indexes <- floor(seq(1, nrow(filtered.mapped.sampled.points), length.out = min(n.samples, nrow(filtered.mapped.sampled.points))))
  real.sampled.points <- filtered.mapped.sampled.points[sample.indexes, ]

  chain.list[[i]] <- coda::as.mcmc(sampled.points[dimensions])
}, mc.cores = min(num.chains, n_cores))

coda.chain.lists <- coda::mcmc.list(chain.list)
#save(coda.chain.lists, file = paste0("~/data/chains/c",http://127.0.0.1:16565/graphics/22c361e1-1c21-482b-8d3b-124b99d51f4c.png
#                                    length(coda.chain.lists),
#                                    "_", length(dimensions),
#                                    "d.RData"))

# METRICS
plot(coda.chain.lists)
traceplot(coda.chain.lists)
geweke.plot(coda.chain.lists)
gelman.plot(coda.chain.lists)
autocorr.plot(coda.chain.lists, lag.max = 50)




