# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)
library(FNN)
library(coda)

result.list <- list()
savedir <- "~/data/chains/"
filename <- "precomputed_chains_50k_env"
datadir <- "/home/mathis/Desktop/semesterarbeit10/data"
result.list[["num.chains"]] <- num.chains <- 5
result.list[["dimensions.list"]] <- dimensions.list <- list(c("PC1", "PC2"),
                                                            c("PC1", "PC2", "PC3"),
                                                            c("PC1", "PC2", "PC3", "PC4"),
                                                            c("PC1", "PC2", "PC3", "PC4","PC5"))
n.samples.per.chain <- 50000
result.list[["burnIn"]] <- burnIn <- FALSE
result.list[["seednumber"]] <- seednumber <- 42



core.distribution <- c(4,4) #max number to be used in each of the two loops

#Needed to stop plotting
plot <- FALSE

# load data
env.data.raster <- USE.MCMC::Worldclim_tmp %>%
  terra::rast( type="xyz")

env.data.sf <- env.data.raster %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))

# fixing the
set.seed(seednumber)

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

# Generate virtual species
virtual.presence.data <- getVirtualSpeciesPresencePoints(env.data = env.data.raster, n.samples = 300)
virtual.presence.points <- virtual.presence.data$sample.points
virtual.presence.points.pc <- terra::extract(rpc$PCs, virtual.presence.points, bind = TRUE) %>%
  sf::st_as_sf()

results.computation <- list()
results.computation <- mclapply(dimensions.list, function(dimensions) {

  # clean data
  env.data.cleaned <- sf::st_drop_geometry(env.with.pc.fs[dimensions])
  env.data.cleaned.subsampled <- sf::st_drop_geometry(env.with.pc.fs.subsampled[dimensions])

  # environment model
  environmental.data.model <- mclust::densityMclust( env.data.cleaned.subsampled, plot = plot)
  summary(environmental.data.model)
  environmental.densities <- mclust::predict.densityMclust(environmental.data.model,  env.data.cleaned.subsampled)
  environment.threshold <- stats::quantile(environmental.densities, 0.01)

  # sample species model
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
  # sample points
  chain.list <- list()
  chain.list <- mclapply(1:num.chains, function(i) {
    sampled.points <- mcmcSampling(dataset = env.with.pc.fs.subsampled,
                                   dimensions = dimensions,
                                   n.sample.points = n.samples.per.chain,
                                   proposalFunction = proposalFunction,
                                   densityFunction = densityFunction,
                                   burnIn = burnIn)
    # mapped.sampled.point.locations <- FNN::get.knnx(env.data.cleaned[dimensions], sampled.points[dimensions],k = 1)
    # mapped.sampled.points <- env.with.pc.fs[mapped.sampled.point.locations$nn.index,]
    # mapped.sampled.points$density <- sampled.points$density
    # mapped.sampled.points$distance <- mapped.sampled.point.locations$nn.dist
    #
    # distance.threshold <- stats::quantile(mapped.sampled.points$distance, 0.95)
    # filtered.mapped.sampled.points <- mapped.sampled.points[mapped.sampled.points$distance < distance.threshold, ]
    #
    #
    # n.samples <- 500
    # sample.indexes <- floor(seq(1, nrow(filtered.mapped.sampled.points), length.out = min(n.samples, nrow(filtered.mapped.sampled.points))))
    # real.sampled.points <- filtered.mapped.sampled.points[sample.indexes, ]

    chain.list[[i]] <- coda::as.mcmc(sampled.points[dimensions])
  }, mc.cores = min(num.chains, core.distribution[1]))

  name <- paste0("mcmc.list.", length(dimensions), "d", num.chains, "c")
  result <- list()
  result[[name]] <- coda::mcmc.list(chain.list)
  return(result)
}, mc.cores = min(num.chains, core.distribution[2]))

result.list[["chains"]] <- unlist(results.computation, recursive = FALSE)
save(result.list, file = paste0(savedir, filename, ".RData"))
