# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)


#Nedded for plotting
par(mfrow = c(1, 1))

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
env.with.pc.fs <- env.with.pc.fs[runif(nrow(env.with.pc.fs)/1, 1, nrow(env.with.pc.fs)),]

#specify the dimension that should be included in the following analysys
dimensions <- c("PC1", "PC2") #, "PC3", "PC4","PC5"

# clean data
env.data.cleaned <- sf::st_drop_geometry(env.with.pc.fs[dimensions])


# environment model
environmental.data.model <- mclust::densityMclust(env.data.cleaned, plot = TRUE)
summary(environmental.data.model)
#env.with.pc.fs$density <- as.integer(mclust::predict.densityMclust(environmental.data.model, env.data.cleaned) < environment.threshold)
env.with.pc.fs$density <- mclust::predict.densityMclust(environmental.data.model, env.data.cleaned)
environment.threshold <- stats::quantile(env.with.pc.fs$density, 0.01)

# sample species model
virtual.presence.data <- getVirtualSpeciesPresencePoints(env.data = env.data.raster, n.samples = 1000)
virtual.presence.points <- virtual.presence.data$sample.points
virtual.presence.points.pc <- terra::extract(rpc$PCs, virtual.presence.points, bind = TRUE) %>%
  sf::st_as_sf()

species.model = mclust::densityMclust(sf::st_drop_geometry(virtual.presence.points.pc[dimensions]), plot = TRUE)
env.with.pc.fs$density.species <- mclust::predict.densityMclust(species.model, env.data.cleaned)
summary(species.model)

#density Function
densityFunction <- mclustDensityFunction(env.model = environmental.data.model, presence.model = species.model,
                                         dim = dimensions, threshold = environment.threshold)

# # set sampling parameters
covariance.proposal.function <-0.5
proposalFunction <- addHighDimGaussian(cov.mat =covariance.proposal.function * diag(length(dimensions)), dim = length(dimensions))


ggplot(data = env.with.pc.fs) +
  geom_sf(aes(color = density)) +  # use density for color mapping
  scale_color_viridis_c() +        # optional: use a nice continuous color scale
  theme_minimal() +                # cleaner theme
  labs(title = "Denisty of the environmental model",
       color = "Density")

ggplot(data = env.with.pc.fs) +
  geom_sf(aes(color = density.species)) +  # use density for color mapping
  scale_color_viridis_c() +        # optional: use a nice continuous color scale
  theme_minimal() +                # cleaner theme
  labs(title = "Denisty of the species model",
       color = "Density")

