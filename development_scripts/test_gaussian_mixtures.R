# Load the mclust package
library(mclust)

# load data

dim <- c("PC1", "PC2")
envData <- USE.MCMC::Worldclim_tmp %>%
  terra::rast( type="xyz")

sfEnvData <- envData %>%
  as.data.frame(xy = TRUE) %>%
  sf::st_as_sf(coords = c("x", "y"))


# Generate the environmental space using PCA
rpc <- rastPCA(envData,  stand = TRUE)


envWithPc <- rpc$PCs %>%
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_join(sfEnvData)

envWithPc <- envWithPc[runif(nrow(envWithPc)/10, 1, nrow(envWithPc)),]

data <- envWithPc[dim] %>%
  st_drop_geometry()


# Now, apply the GMM using mclust
gmm_model <- densityMclust(data)

# View the model summary
summary(gmm_model)

# Plot the clustering result
plot(gmm_model)

# Display the best number of clusters and model type
cat("Best number of clusters: ", gmm_model$G, "\n")
cat("Best model type: ", gmm_model$modelName, "\n")

