library(gganimate)
library(tcltk)
library(coda)

animateChain <- function(chain) {
  p <- ggplot(chain, aes(x = PC1, y = PC2, color = PC3)) +
    geom_path() +
    labs(title = 'Chain runntime',
         subtitle = 'Point number: {frame_time}')

  anim <- p + transition_time(year) + ease_aes('linear') +
    #transition_reveal(year) +
    exit_shrink() +
    shadow_trail(alpha = 0.3)

}

animateChain2d <- function(chain) {
  p <- ggplot(chain, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.7) +
    # geom_line() +
    labs(title = 'Chain runntime',
         subtitle = 'Point number: {frame_time}')

  anim <- p + transition_time(year) + ease_aes('linear') +
    exit_shrink() +
    shadow_trail(alpha = 0.3)

}
fps <- 24
duration <- 10
n.frames <- fps * duration

load("/home/mathis/Desktop/semesterarbeit10/data/chains/precomputed_chains_50k_env.RData")
load("/home/mathis/Desktop/semesterarbeit10/data/chains/c4_3d.RData")
df <- sampled.points
df <- as.data.frame(result.list$chains$mcmc.list.2d5c[[1]])[1:20000,]
#df <- as.data.frame(coda.chain.lists[[1]])[1:10000,]
df$year <- floor(1:nrow(df) / n.frames)
anim <- animateChain2d(df)
animate(anim, renderer = gifski_renderer(), width = 800, height = 800, fps = fps, duration = duration)

anim_save("../plots/gifs/my_animation.gif", animation = gganimate::last_animation())
