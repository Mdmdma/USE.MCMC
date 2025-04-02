# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)

# load data
envData <- USE.MCMC::Worldclim_tmp %>%
  terra::rast( type="xyz")

sfEnvData <- envData %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))

# Create virtual species
set.seed(123)
random.sp <- virtualspecies::generateRandomSp(envData,
                                              convert.to.PA = FALSE,
                                              species.type = "additive",
                                              realistic.sp = TRUE,
                                              plot = FALSE)

# Reclassify suitability raster using a probability conversion rule
new.pres <- virtualspecies::convertToPA(x=random.sp,
                                        beta=0.55,
                                        alpha = -0.05, plot = FALSE)

# Sample true occurrences
presence.points <- virtualspecies::sampleOccurrences(new.pres,
                                                     n = 300, # The number of points to sample
                                                     type = "presence-absence",
                                                     sample.prevalence = 0.6,
                                                     detection.probability = 1,
                                                     correct.by.suitability = TRUE,
                                                     plot = TRUE)

# Generate a presence-only data set
myPres <- presence.points$sample.points[which(presence.points$sample.points$Observed==1), c("x", "y", "Observed")]
myPres <- st_as_sf(myPres, coords=c("x", "y"), crs=4326)

# Generate the environmental space using PCA
rpc <- rastPCA(envData,  stand = TRUE)


envWithPc <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(sfEnvData)

# subsample env space to speed up the proces
# envWithPc <- envWithPc[runif(nrow(envWithPc)/10, 1, nrow(envWithPc)),]

# sample points
sampled.points <- mcmcSampling(dataset = envWithPc, dimensions = c("PC1", "PC2", "PC3"), n.sample.points = 1000)
par(mfrow = c(2, 1))
plot_points_with_lines(sampled.points, c("PC1", "PC2"))
plot(envWithPc$PC1, envWithPc$PC2)
par(mfrow = c(2, 2))
hist(envWithPc$PC1)
hist(envWithPc$PC2)
hist(sampled.points$PC1)
hist(sampled.points$PC2)




