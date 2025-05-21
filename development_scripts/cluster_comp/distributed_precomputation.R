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
datadir <- "/home/mathis/Desktop/semesterarbeit10/"
result.list[["num.chains"]] <- num.chains <- 10
result.list[["dimensions.list"]] <- dimensions.list <- list(c("PC1", "PC2"),
                                                            c("PC1", "PC2", "PC3"),
                                                            c("PC1", "PC2", "PC3", "PC4"),
                                                            c("PC1", "PC2", "PC3", "PC4","PC5"))

result.list[["dimensions.list"]] <- dimensions.list <- list(c("wc2.1_10m_bio_3", "wc2.1_10m_bio_4", "wc2.1_10m_bio_9", "wc2.1_10m_bio_14", "wc2.1_10m_bio_15"))

n.samples.per.chain <- 50000
result.list[["burnIn"]] <- burnIn <- TRUE
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
env.data.raster <- c(env.data.raster, rpc$PCs) # TODO change nameing to be better and not need this weired positioning
virtual.presence.points.pc <- terra::extract(env.data.raster, virtual.presence.points, bind = TRUE) %>%
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
  species.densities <- species.model$density
  species.cutoff.threshold <- stats::quantile(species.densities, 0.9)

  #density Function
  densityFunction <- mclustDensityFunction(env.model = environmental.data.model,
                                           species.model = species.model,
                                           dim = dimensions,
                                           threshold = environment.threshold,
                                           species.cutoff.threshold = species.cutoff.threshold)

  # # set sampling parameters
  covariance.scaling <-0.075
  covariance.matrix <- stats::cov(sf::st_drop_geometry(env.with.pc.fs)[dimensions])
  proposalFunction <- addHighDimGaussian(cov.mat =covariance.scaling * covariance.matrix,
                                         dim = length(dimensions))
  # sample points
  chain.list <- list()
  chain.list <- mclapply(1:num.chains, function(i) {
    capture.output({
      sampled.points <- mcmcSampling(dataset = env.with.pc.fs.subsampled,
                                     dimensions = dimensions,
                                     n.sample.points = n.samples.per.chain,
                                     proposalFunction = proposalFunction,
                                     densityFunction = densityFunction,
                                     burnIn = burnIn,
                                     verbose = FALSE)
    })
    chain.list[[i]] <- coda::as.mcmc(sampled.points[dimensions])
  }, mc.cores = min(num.chains, core.distribution[1]))

  name <- paste0("mcmc.list.", length(dimensions), "d", num.chains, "c")
  result <- list()
  result[[name]] <- coda::mcmc.list(chain.list)
  return(result)
}, mc.cores = min(num.chains, core.distribution[2]))

result.list[["chains"]] <- unlist(results.computation, recursive = FALSE)
save(result.list, file = paste0(savedir, filename, ".RData"))
