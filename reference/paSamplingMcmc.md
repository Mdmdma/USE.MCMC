# paSamplingMcmc is a near drop in replacement for paSampling from the original USE package, that allows to perform a Gaussian mixture based pseudo absence sampling using a markov. In a first step a density function is constructed using a GMM fitted to the environment as a limit to the sampling space and a GMM fitted on the target species as a way to evade regions associated with the presence.

paSamplingMcmc is a near drop in replacement for paSampling from the
original USE package, that allows to perform a Gaussian mixture based
pseudo absence sampling using a markov. In a first step a density
function is constructed using a GMM fitted to the environment as a limit
to the sampling space and a GMM fitted on the target species as a way to
evade regions associated with the presence.

## Usage

``` r
paSamplingMcmc(
  env.data.raster = NULL,
  pres = NULL,
  n.samples = 300,
  chain.length = 10000,
  verbose = FALSE,
  dimensions = c("PC1", "PC2"),
  burnIn = 1000,
  covariance.correction = 1,
  precomputed.pca = NULL,
  seed.number = 42,
  n.neighbors.for.statistics = 2,
  low.end.of.inclueded.points = 100,
  high.end.of.included.points = 5,
  environmental.cutof.percentile = 0.001,
  species.cutoff.threshold = 0.95,
  plot_proc = FALSE,
  num.chains = 1,
  num.cores = 1
)
```

## Arguments

- env.data.raster:

  Terra raster containing the environment

- pres:

  Sf dataframe containing the presence locations

- n.samples:

  number of samples that should be put out

- chain.length:

  number of points that are sampled for the chain

- verbose:

  If true the function gives updates on the current state of the chain

- dimensions:

  vector containg the names of the dimensions that should be included

- burnIn:

  Integer, sets the number of steps per adaptive burnin cycle. If 0 the
  burnin is skipped

- covariance.correction:

  Integer, sets the inital value of the covariance correction

- precomputed.pca:

  If rastPCA has already been evoked, it the result of it can be passed
  here to not recompute

- seed.number:

  seednumber used to get repeatable results

- n.neighbors.for.statistics:

  number of neighbors used to calculate the maximal sensible distance to
  real points that should be included

- low.end.of.inclueded.points:

  Sets the range of points included in the threshold computation

- high.end.of.included.points:

  Sets the range of points included in the threshold computation

- environmental.cutof.percentile:

  sets the percentile of the environment GMM that is excluded from the
  space that can be visited by the chain

- species.cutoff.threshold:

  sets the percentile of the species presence GMM that is included in
  the space that can be visited by the chain

- plot_proc:

  If true the function returns plots the progress

- num.chains:

  Number of chains from which samples should be picked

- num.cores:

  Number of cores available for parallelization of the multi-chain
  computation

## Value

dataframe containing the sampled points
