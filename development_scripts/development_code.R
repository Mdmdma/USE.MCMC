# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)

# Define data directory and download WorldClim data
dataDir <- "/home/mathis/Desktop/semesterarbeit10/data/demo"
Worldclim <- geodata::worldclim_global(var = 'bio', res = 2.5, path = dataDir)
names(Worldclim) <- paste0("bio", 1:length(names(Worldclim)))

# Load environmental data and crop to European extent
# In the actual script you would use:
# Worldclim <- geodata::worldclim_global(var='bio', res=10, path=getwd())
envData <- terra::crop(Worldclim, terra::ext(-12, 25, 36, 60))

# In the vignette, they're using pre-loaded data:
# envData <- USE.MCMC::Worldclim_tmp
# envData <- terra::rast(envData, type="xy")

# Select random variables
# myRandNum <- sample(1:19, size=5, replace = FALSE)
# envData <- envData[[myRandNum]]

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
rpc <- rastPCA(envData, stand = TRUE)
dt <- na.omit(as.data.frame(rpc$PCs[[c("PC1", "PC2")]], xy = TRUE))
dt <- sf::st_as_sf(dt, coords = c("PC1", "PC2"))

# Find optimal resolution for sampling grid
# myRes <- USE.MCMC::optimRes(sdf=dt,
#                     grid.res=c(1:10),
#                     perc.thr = 20,
#                     showOpt = TRUE,
#                     cr=5)

# In the vignette, they're using a pre-computed result:
myRes <- list()
myRes$Opt_res <- 5
myRes$Opt_res  # Display optimal resolution

# Uniform sampling of the environmental space
myObs <- USE.MCMC::uniformSampling(sdf=dt,
                                   grid.res=myRes$Opt_res,
                                   n.tr = 5,
                                   sub.ts = TRUE,
                                   n.ts = 2,
                                   plot_proc = FALSE)

# View sampled observations
head(myObs$obs.tr)

# Visualize the coordinates (PC-scores) of sampled observations
env_pca <- c(rpc$PCs$PC1, rpc$PCs$PC2)
env_pca <- na.omit(as.data.frame(env_pca))

# Plot density of PC1 scores
ggplot(env_pca, aes(x=PC1))+
  geom_density(aes(color="Environment"), linewidth=1 )+
  geom_density(data=data.frame(st_coordinates(myObs$obs.tr)),
               aes(x=X,  color="Uniform"), linewidth=1)+
  scale_color_manual(name=NULL,
                     values=c('Environment'='#1E88E5', 'Uniform'='#D81B60'))+
  labs(y="Density of PC-scores")+
  ylim(0,1)+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=14),
        legend.text=element_text(size=12))

# Plot density of PC2 scores
ggplot(env_pca, aes(x=PC2))+
  geom_density(aes(color="Environment"), linewidth=1 )+
  geom_density(data=data.frame(st_coordinates(myObs$obs.tr)),
               aes(x=Y,  color="Uniform"), size=1)+
  scale_color_manual(name=NULL,
                     values=c('Environment'='#1E88E5', 'Uniform'='#D81B60'))+
  labs(y="Density of PC-scores")+
  ylim(0,1)+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=14),
        legend.text=element_text(size=12))

# Uniform sampling of pseudo-absences
myGrid.psAbs <- USE.MCMC::paSampling(env.rast=envData,
                                     pres=myPres,
                                     thres=0.75,
                                     H=NULL,
                                     grid.res=as.numeric(myRes$Opt_res),
                                     n.tr = 5,
                                     prev=0.3,
                                     sub.ts=TRUE,
                                     n.ts=5,
                                     plot_proc=FALSE,
                                     verbose=FALSE)

# Visualize PC1 density of pseudo-absences
ggplot(env_pca, aes(x=PC1))+
  geom_density(aes(color="Environment"), linewidth=1 )+
  geom_density(data=data.frame(st_coordinates(myGrid.psAbs$obs.tr)),
               aes(x=X,  color="Uniform"), linewidth=1)+
  geom_density(data=terra::extract(c(rpc$PCs$PC1, rpc$PCs$PC2), myPres, df=TRUE),
               aes(x=PC1, color="Presence"), linewidth=1 )+
  scale_color_manual(name=NULL,
                     values=c('Environment'='#1E88E5', 'Uniform'='#D81B60', "Presence"="black"))+
  labs(y="Density of PC-scores")+
  ylim(0,1)+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=14),
        legend.text=element_text(size=12))

# Visualize PC2 density of pseudo-absences
ggplot(env_pca, aes(x=PC2))+
  geom_density(aes(color="Environment"), linewidth=1 )+
  geom_density(data=data.frame(st_coordinates(myGrid.psAbs$obs.tr)),
               aes(x=Y,  color="Uniform"), linewidth=1)+
  geom_density(data=terra::extract(c(rpc$PCs$PC1, rpc$PCs$PC2), myPres, df=TRUE),
               aes(x=PC2, color="Presence"), linewidth=1 )+
  scale_color_manual(name=NULL,
                     values=c('Environment'='#1E88E5', 'Uniform'='#D81B60', "Presence"="black"))+
  labs(y="Density of PC-scores")+
  ylim(0,1)+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=14),
        legend.text=element_text(size=12))

# Visualize geographic coordinates of sampled pseudo-absences
ggplot()+
  tidyterra::geom_spatraster(data = new.pres$pa.raster)+
  scale_fill_viridis_c(na.value = "transparent", breaks=c(0,1)) +
  geom_sf(data=myPres,
          aes(color= "Presences"),
          alpha=1, size=2, shape= 19)+
  geom_sf(data=st_as_sf(st_drop_geometry(myGrid.psAbs$obs.tr),
                        coords = c("x","y"), crs=4326),
          aes(color="Pseudo-absences"),
          alpha=0.8, size=2, shape = 19 )+
  scale_colour_manual(name=NULL,
                      values=c('Presences'='steelblue','Pseudo-absences'='#A41616'))+
  labs(x="Longitude",
       y="Latitude",
       fill="Virtual species")+
  theme_light()+
  theme(legend.position = "bottom",
        legend.background=element_blank(),
        legend.box="vertical",
        panel.grid = element_blank(),
        text = element_text(size=14),
        legend.text=element_text(size=14),
        aspect.ratio = 1,
        panel.spacing.y = unit(2, "lines"))

# Inspect effect of different threshold values on the environmental sub-space
USE.MCMC::thresh.inspect(env.rast=envData,
                         pres=myPres,
                         thres=c(0.1, 0.25, 0.5, 0.75, 0.9),
                         H=NULL
)

