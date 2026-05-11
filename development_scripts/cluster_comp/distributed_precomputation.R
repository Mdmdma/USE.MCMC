# Distributed precomputation of MCMC chains across PC-dimension counts.
#
# Mirrors the setup that paSamplingMcmc does internally, but stops at the raw
# coda::mcmc.list level so the chains can be inspected with convergence
# diagnostics (Gelman-Rubin, autocorrelation, ...). The setup pipeline routes
# through the same optimized paths that the library now uses (terra-based
# combined raster, gated verbose output, partial sort in distance threshold).

library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(parallel)
library(FNN)
library(coda)
library(mclust)
library(optparse)

`%>%` <- magrittr::`%>%`

option_list <- list(
  make_option(c("-s", "--seed"), type = "integer", help = "seed number", default = 42),
  make_option(c("--pca"),         type = "logical", help = "if false sample in original bioclim variables", default = TRUE),
  make_option(c("--burnin"),      type = "integer", help = "burn-in iterations per adaptive cycle",         default = 1000L),
  make_option(c("--chain-length"),type = "integer", help = "samples per chain",                              default = 50000L),
  make_option(c("--num-chains"),  type = "integer", help = "chains per dimension count",                     default = 10L),
  make_option(c("--env-quantile"),type = "double",  help = "environment density cutoff quantile",            default = 0.001),
  make_option(c("--species-quantile"), type = "double", help = "species density cutoff quantile",            default = 0.95)
)
opt <- parse_args(OptionParser(option_list = option_list))

savedir <- Sys.getenv("USE_MCMC_CHAINS_DIR", unset = path.expand("~/data/chains/"))
if (!dir.exists(savedir)) dir.create(savedir, recursive = TRUE)

result.list <- list()
result.list[["num.chains"]] <- num.chains <- opt$`num-chains`
result.list[["burnIn"]]     <- burnIn     <- as.integer(opt$burnin)
result.list[["seed.number"]] <- seed.number <- opt$seed
n.samples.per.chain <- as.integer(opt$`chain-length`)

if (opt$pca) {
  result.list[["dimensions.list"]] <- dimensions.list <- list(
    c("PC1", "PC2"),
    c("PC1", "PC2", "PC3"),
    c("PC1", "PC2", "PC3", "PC4"),
    c("PC1", "PC2", "PC3", "PC4", "PC5")
  )
} else {
  result.list[["dimensions.list"]] <- dimensions.list <- list(
    c("wc2.1_10m_bio_3", "wc2.1_10m_bio_4", "wc2.1_10m_bio_9",
      "wc2.1_10m_bio_14", "wc2.1_10m_bio_15")
  )
}

file.prefix <- "precomputed_chains_n"
filename    <- paste0(file.prefix, n.samples.per.chain, "_c", num.chains, "_pca", as.character(opt$pca))
core.distribution <- c(4, 4) # max cores per (inner-chains, outer-dimensions) loops

# load data and ensure a SpatRaster (matches the new paSamplingMcmc setup)
env.data.raster <- USE.MCMC::Worldclim_tmp %>%
  terra::rast(type = "xyz")

set.seed(seed.number)
rpc <- rastPCA(env.data.raster, stand = TRUE)

# Combined env+PC dataframe via the same fast path the library uses.
env.with.pc.sf <- terra::as.data.frame(c(env.data.raster, rpc$PCs), xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y"))

# subsample env space to keep mclust density fitting tractable
env.with.pc.sf.subsampled <- env.with.pc.sf[
  stats::runif(min(nrow(env.with.pc.sf), 2000), 1, nrow(env.with.pc.sf)), ]

# Generate virtual species and attach PC coords once
virtual.presence.data   <- getVirtualSpeciesPresencePoints(env.data = env.data.raster, n.samples = 300)
virtual.presence.points <- virtual.presence.data$sample.points
env.combined.raster     <- c(env.data.raster, rpc$PCs)
virtual.presence.points.pc <- terra::extract(env.combined.raster, virtual.presence.points, bind = TRUE) %>%
  sf::st_as_sf()

results.computation <- mclapply(dimensions.list, function(dimensions) {
  env.data.cleaned            <- sf::st_drop_geometry(env.with.pc.sf[dimensions])
  env.data.cleaned.subsampled <- sf::st_drop_geometry(env.with.pc.sf.subsampled[dimensions])

  environmental.data.model <- mclust::densityMclust(env.data.cleaned.subsampled,
                                                    plot = FALSE, verbose = FALSE)
  environmental.densities <- mclust::predict.densityMclust(environmental.data.model,
                                                           env.data.cleaned)
  environment.threshold <- stats::quantile(environmental.densities, opt$`env-quantile`)

  species.model <- mclust::densityMclust(
    sf::st_drop_geometry(virtual.presence.points.pc[dimensions]),
    plot = FALSE, verbose = FALSE)
  species.densities <- species.model$density
  species.cutoff <- stats::quantile(species.densities, opt$`species-quantile`)

  densityFunction <- mclustDensityFunction(env.model = environmental.data.model,
                                           species.model = species.model,
                                           dim = dimensions,
                                           threshold = environment.threshold,
                                           species.cutoff.threshold = species.cutoff)

  covariance.matrix <- stats::cov(sf::st_drop_geometry(env.with.pc.sf)[dimensions])
  proposalFunction  <- addHighDimGaussian(cov.mat = 0.075 * covariance.matrix,
                                          dim = length(dimensions))

  chain.list <- mclapply(seq_len(num.chains), function(i) {
    sampled.points <- mcmcSampling(dataset = env.with.pc.sf.subsampled,
                                   dimensions = dimensions,
                                   n.sample.points = n.samples.per.chain,
                                   proposalFunction = proposalFunction,
                                   densityFunction = densityFunction,
                                   burnIn = burnIn,
                                   verbose = FALSE)
    coda::as.mcmc(sampled.points[dimensions])
  }, mc.cores = min(num.chains, core.distribution[1]))

  result.name <- paste0("mcmc.list.", length(dimensions), "d", num.chains, "c")
  out <- list()
  out[[result.name]] <- coda::mcmc.list(chain.list)
  out
}, mc.cores = min(length(dimensions.list), core.distribution[2]))

result.list[["chains"]] <- unlist(results.computation, recursive = FALSE)
save(result.list, file = file.path(savedir, paste0(filename, ".RData")))
