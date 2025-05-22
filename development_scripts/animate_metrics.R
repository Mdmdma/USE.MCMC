library(coda)
library(ggplot2)
library(gganimate)
library(tcltk)
library(tidyr)

create.frame.indexes <- function(df, n.frames){
  total_rows <- nrow(df)
  rows_per_frame <- floor(seq(1, total_rows, length.out = n.frames + 1))

  result_rows <- sum(diff(rows_per_frame))
  result <- data.frame(matrix(NA, nrow = result_rows, ncol = ncol(df) + 1))
  names(result) <- c(names(df), "frame.number")

  current_row <- 1
  for (i in 1:n.frames) {
    # Calculate start and end rows for this frame
    start_row <- 1  # Always start from the beginning
    end_row <- rows_per_frame[i + 1]

    # Get rows for this frame
    frame_size <- end_row - start_row + 1
    result_indices <- current_row:(current_row + frame_size - 1)

    # Copy data and set frame number
    result[result_indices, 1:ncol(df)] <- df[start_row:end_row, ]
    result[result_indices, "frame.number"] <- i

    # Update current position in result dataframe
    current_row <- current_row + frame_size
  }

  return(result)
}


filepath <- "/home/mathis/Desktop/semesterarbeit10/data/chains/precomputed_chains_50k_pca_s42.RData"
load(filepath)
fps <- 10
duration <- 1
n.frames <- fps * duration

variable <- "PC1"
chainset <- result.list$chains[[1]]
df.mc <- as.data.frame(chainset[[1]])[variable]
names(df.mc) <- "chain 1"


for (chain.number in 2:length(chainset)) {
  chain.name <- paste0("chain ", chain.number)
  print(chain.name)
  df.mc[chain.name] <- as.data.frame(chainset[[chain.number]])[variable]
  df.mc <- as.data.frame(df.mc)
}

frames <- create.frame.indexes(df.mc, n.frames = n.frames)
frames.long <- pivot_longer(frames,
                          cols = names(df.mc),
                          names_to = "chain",
                          values_to = "value")

plot <- ggplot(frames.long, aes(x = value, color = chain)) +
          geom_density() +
          labs(title = "Density plot",
               subtitle = 'Chain length: {frame_time}')

anim <- plot + transition_manual(frame.number)

animated <- animate(anim, renderer = gifski_renderer(), width = 600, height = 400, fps = fps, duration = duration)

anim_save("../plots/gifs/density2d_pc1_10c_1s.gif", animation = animated)








