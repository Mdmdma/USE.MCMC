library(entropy)

plot_entropy <- function(chain.list){
  interval <-  floor(seq(1, 50000, length.out = 100))
  results <- data.frame(interval = numeric(0), Hn = numeric(0), chain_id = integer(0))
  for (i in seq_along(chain.list)) {
    chain <- chain.list[i]

    # Calculate entropy for each interval and store it in a temporary data frame
    H <- sapply(interval, function(nrow) entropy(as.matrix(chain)[1:nrow, ]))

    # Add the results for the current chain to the results data frame
    results <- rbind(results, data.frame(interval = interval, H = H, chain_id = rep(i, length(interval))))
  }

  # Now create the plot
  plot <- ggplot(results, aes(x = interval, y = H, color = as.factor(chain_id))) +
    geom_line() +
    labs(color = "Chain",
         y = "Entopy",
         x = "Points included until",
         title = paste0("Entropy along the Chain in ",dim(as.matrix(chain.list[1]))[2],"d"))+
    theme_minimal() +
    ylim(0,max(100, min(results)))
  plot
}


plot_KLdiff <- function(chain.list){
  interval <-  floor(seq(1, 50000, length.out = 100))
  results <- data.frame(interval = numeric(0), Hn = numeric(0), chain_id = integer(0))
  for (i in seq_along(chain.list)) {
    chain <- chain.list[i]

    # Calculate entropy for each interval and store it in a temporary data frame
    H <- sapply(interval, function(nrow) FNN::KL.divergence(as.matrix(chain)[1:nrow, ], as.matrix(chain))[3])

    # Add the results for the current chain to the results data frame
    results <- rbind(results, data.frame(interval = interval, H = H, chain_id = rep(i, length(interval))))
  }

  # Now create the plot
  plot <- ggplot(results, aes(x = interval, y = H, color = as.factor(chain_id))) +
    geom_line() +
    labs(color = "Chain",
         y = "KL divergence",
         x = "Points included until",
         title = paste0("KL divergence along the chain in ",dim(as.matrix(chain.list[1]))[2],"d"))+
    theme_minimal() +
    ylim(0,max(100, min(results)))
  plot
}


#load precomputed chains from files
env5d <- new.env()
env4d <- new.env()
env3d <- new.env()
env2d <- new.env()

# Load each file into its own environment
load("/run/user/1000/gvfs/sftp:host=euler.ethz.ch/cluster/home/merler/data/chains/c4_5d.RData", envir = env5d)
load("/run/user/1000/gvfs/sftp:host=euler.ethz.ch/cluster/home/merler/data/chains/c4_4d.RData", envir = env4d)
load("/run/user/1000/gvfs/sftp:host=euler.ethz.ch/cluster/home/merler/data/chains/c4_3d.RData", envir = env3d)
load("/run/user/1000/gvfs/sftp:host=euler.ethz.ch/cluster/home/merler/data/chains/c4_2d.RData", envir = env2d)

d5 <- env5d$coda.chain.lists
d4 <- env4d$coda.chain.lists
d3 <- env3d$coda.chain.lists
d2 <- env2d$coda.chain.lists

gelman.plot(d2)
gelman.plot(d3)
gelman.plot(d4)
gelman.plot(d5)

chain.list <- d5
plot(chain.list)
traceplot(d2)
geweke.plot(d5)
gelman.plot(d2)
autocorr.plot(d2, lag.max = 10000)

plot_entropy(d2)

if (TRUE){
  #analyse predictive model performance
  chain <- d3[[1]]
  species.raster <- virtual.presence.data$original.distribution.raster
  vs.model.df <- as.data.frame(lapply(species.raster, as.numeric))
  random_expected <- mean(vs.model.df$lyr.1)

  sampled.points <- as.data.frame(chain)
  mapped.sampled.point.locations <- FNN::get.knnx(env.data.cleaned[dimensions], sampled.points[dimensions],k = 1)
  mapped.sampled.points <- env.with.pc.fs[mapped.sampled.point.locations$nn.index,]
  mapped.sampled.points$density <- sampled.points$density
  mapped.sampled.points$distance <- mapped.sampled.point.locations$nn.dist

  distance.threshold <- stats::quantile(mapped.sampled.points$distance, 0.95)
  filtered.mapped.sampled.points <- mapped.sampled.points[mapped.sampled.points$distance < distance.threshold, ]


  n.samples <- 10000
  sample.indexes <- floor(seq(1, nrow(filtered.mapped.sampled.points), length.out = min(n.samples, nrow(filtered.mapped.sampled.points))))
  real.sampled.points <- filtered.mapped.sampled.points[sample.indexes, ]
  # TODO check where the coordinate system gets lost
  st_crs(real.sampled.points) <- 4326

  vs.model.prediction.on.sample.points <- terra::extract(species.raster, real.sampled.points)
  model_perf <- mean(vs.model.prediction.on.sample.points$lyr.1)
  print(model_perf)
}

