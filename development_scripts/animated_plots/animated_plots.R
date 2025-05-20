library(ggplot2)
library(gganimate)
library(gapminder)
library(tidyr)

# If your data is in wide format, reshape it first
# Replace parameter_names with your actual parameter column names
df <- sampled.points
df$iteration <- floor(1:nrow(df) / 100)
mcmc_long <- pivot_longer(df,
                          cols = c("PC1", "iteration"),
                          names_to = "parameter",
                          values_to = "value")

# Create the animation
p <- ggplot(mcmc_long, aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~parameter, scales = "free") +
  labs(title = "Posterior Density Evolution",
       subtitle = "Iteration: {frame_time}") +
  ylim(-5,5) +
  theme_minimal()

# Add animation
anim <- p +
  transition_time(iteration) +
  ease_aes('linear')

# Render with appropriate settings
gganimate::animate(anim, nframes = 100, fps = 10, width = 800, height = 500)
