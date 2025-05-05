library(mclust)
library(tidyterra)
# Gaussian mixtures basics
set.seed(123)
par(mfrow = c(1, 3))
# Parameters
n <- 1000               # total number of data points
pi <- c(1, 1, 1)        # mixing proportions
mu <- c(-2, 3, 1)           # means of the Gaussians
sigma <- c(1, 0.5, 2)       # standard deviations

# Simulate component membership
z <- sample(1:3, size = n, replace = TRUE, prob = pi)

# Generate data
x <- rnorm(n, mean = mu[z], sd = sigma[z])

plot(x, main = "Datapoints")
# Visualize the histogram
#plot(density(x),  main = "Data density")
hist(x, main = "Data histogram")

model <- Mclust(x, plot = FALSE)  # G = number of components

summary(model)

# Plot the fitted model
plot(model, what = "density", data = x)

# Plot histogram of data
#hist(x, breaks = 50, probability = TRUE, col = "lightgray", main = "Gaussian Mixture Model with Components")


# Define x-range for plotting densities
x_vals <- seq(min(x), max(x), length.out = 1000)

# Extract estimated parameters
props <- model$parameters$pro
means <- model$parameters$mean
sds <- sqrt(model$parameters$variance$sigmasq)

# Plot each Gaussian component
colors <- c("red", "blue", "green")
for (k in 1:3) {
  lines(x_vals, props[k] * dnorm(x_vals, mean = means[k], sd = sds[k]), col = colors[k], lwd = 2)
}

# Optional: add overall mixture density
mixture_density <- rowSums(sapply(1:3, function(k) props[k] * dnorm(x_vals, mean = means[k], sd = sds[k])))
lines(x_vals, mixture_density, col = "black", lwd = 2, lty = 2)

legend("topright", legend = c("C 1", "C 2", "C 3", "Mix"),
       col = c("red", "blue", "green", "black"), lty = c(1, 1, 1, 2), lwd = 2)
title("Fitted model")


#!! needs to be run after mcmc_test_script.R !!
# Gaussian Mixtures in multiple dimensions


par(mfrow = c(1, 3))
plot(env.data.cleaned)
title("Environmental data")
plotDensityMclust2(environmental.data.model, nlevels = 10, levels = seq(0.001, 0.5, length.out=50),
                   xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                   ylim = c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)))
title("Predicted env density")
plotDensityMclust2(species.model,
                   xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                   ylim = c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)))
title("Predicted species density")

# model density
plotDensityLines(dataset = real.sampled.points[1,],
                 lines = FALSE, cols = c("PC1", "PC2", "PC3"),
                 density = TRUE, densityFunction = densityFunction, resolution = 100,
                 xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                 ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                 title = "Density Function"
                 )
plotDensityLines(dataset = sampled.points[1:200, ],
                 lines = TRUE, cols = c("PC1", "PC2", "PC3"),
                 density = TRUE, densityFunction = densityFunction, resolution = 100, species = virtual.presence.points.pc,
                 xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                 ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                 title = "Density with part of a chain"
                  )
plotDensityLines(dataset = real.sampled.points,
                 lines = FALSE, cols = c("PC1", "PC2", "PC3"),
                 density = TRUE, densityFunction = densityFunction, resolution = 100, species = virtual.presence.points.pc,
                 xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                 ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                 title = "Density with sampled points"
)

raster_layer <- rpc$PCs$PC2
output_file <- "xmy_raster_plot.png"

# Open PNG device with specific dimensions
# Width and height in pixels (adjust as needed)
width_px <- 200
height_px <- 140
png(filename = output_file, width = width_px, height = height_px)

# --- Create the plot using ggplot2 and tidyterra ---
p <- ggplot() +
  # First add the raster data
  geom_spatraster(data = raster_layer) +
  # Then apply the color scale with NA set to white
  scale_fill_viridis_c(
    na.value = "white"
  ) +
  theme_void() + # Use theme_void for a white background and no axis text/ticks
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "pt"), # Set all margins to 0 with pt units
    panel.border = element_blank(), # Remove panel border
    panel.spacing = unit(0, "pt"), # Remove panel spacing
    panel.background = element_rect(fill = "white", colour = NA), # Ensure panel background is white
    plot.background = element_rect(fill = "white", colour = NA), # Ensure plot background is white
    legend.position = "none" # Remove the legend
  ) +
  # Add this to ensure the plot expands to fill the entire space
  coord_sf(expand = FALSE)

# Print the plot to the device
print(p)

# Close the device - this is essential to save the file
dev.off()

#plot model
par(mfrow = c(1, 2))
plotDensityMclust2(environmental.data.model, nlevels = 10, levels = seq(0.001, 0.5, length.out=50),
                   xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                   ylim = c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)))

plotDensityMclust2(species.model,
                   xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                   ylim = c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)))

plotDensityLines(dataset = real.sampled.points[1:200,],
                 lines = TRUE, cols = c("PC1", "PC2", "PC3"),
                 density = FALSE, densityFunction = densityFunction, resolution = 100, species = virtual.presence.points.pc,
                 xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                 ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                 title = NULL
)

plotDensityLines(dataset = real.sampled.points[1:200,],
                 lines = FALSE, cols = c("PC1", "PC2", "PC3"),
                 density = TRUE, densityFunction = densityFunction, resolution = 100, species = virtual.presence.points.pc,
                 xlim = c(min(env.with.pc.fs$PC1), max(env.with.pc.fs$PC1)),
                 ylim =c(min(env.with.pc.fs$PC2), max(env.with.pc.fs$PC2)),
                 title = NULL
)


