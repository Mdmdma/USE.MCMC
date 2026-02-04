# Package index

## All functions

- [`SpatialProba()`](https://mdmdma.github.io/USE.MCMC/reference/SpatialProba.md)
  : Virtual species probability of occurrence
- [`Worldclim_tmp`](https://mdmdma.github.io/USE.MCMC/reference/Worldclim_tmp.md)
  : A subset of WorldClim bioclimatic variables
- [`addHighDimGaussian()`](https://mdmdma.github.io/USE.MCMC/reference/addHighDimGaussian.md)
  : addHighDimGaussian
- [`alwaysOne()`](https://mdmdma.github.io/USE.MCMC/reference/alwaysOne.md)
  : AlwaysOne
- [`getVirtualSpeciesPresencePoints()`](https://mdmdma.github.io/USE.MCMC/reference/getVirtualSpeciesPresencePoints.md)
  : Function to streamline the generation and sampling of a virtual
  species
- [`mapBackOnRealPoints()`](https://mdmdma.github.io/USE.MCMC/reference/mapBackOnRealPoints.md)
  : Map Back on Real Points searches the closest point in the dataset
  regarding the given point and the dimensions given
- [`maxResNn()`](https://mdmdma.github.io/USE.MCMC/reference/maxResNn.md)
  : maxResNN is a function that can be used to compute a reasonable grid
  resolution for nearest neighbor based uniform sampling. The core idea
  behind its working principle is that we want to expect a grid cell to
  contain points if it overlaps with the environment. This
  implementation looks at the low density regions using distance to n
  neighbors as a proxy. The approach assumes that both coordinates have
  similar range, as the axes are not weighted when computing Neighbors
  and converting from distances to number of grid cells If multiple
  points from the same grid cell should be sampled, the number of
  neighbors included in the computation should be set accordingly
- [`mclustDensityFunction()`](https://mdmdma.github.io/USE.MCMC/reference/mclustDensityFunction.md)
  : Helper to create a Density function that uses mclust Gaussian
  mixtures
- [`mcmcSampling()`](https://mdmdma.github.io/USE.MCMC/reference/mcmcSampling.md)
  : MCMC sampling from a given dataset
- [`optimRes()`](https://mdmdma.github.io/USE.MCMC/reference/optimRes.md)
  : Get optimal resolution of the sampling grid
- [`optimalDistanceThresholdNn()`](https://mdmdma.github.io/USE.MCMC/reference/optimalDistanceThresholdNn.md)
  : Function to calculate the optimal distance threshold for nearest
  neighborhood sampling.
- [`paSampling()`](https://mdmdma.github.io/USE.MCMC/reference/paSampling.md)
  : Sampling pseudo-absences for the training and testing datasets.
- [`paSamplingMcmc()`](https://mdmdma.github.io/USE.MCMC/reference/paSamplingMcmc.md)
  : paSamplingMcmc is a near drop in replacement for paSampling from the
  original USE package, that allows to perform a Gaussian mixture based
  pseudo absence sampling using a markov. In a first step a density
  function is constructed using a GMM fitted to the environment as a
  limit to the sampling space and a GMM fitted on the target species as
  a way to evade regions associated with the presence.
- [`paSamplingNn()`](https://mdmdma.github.io/USE.MCMC/reference/paSamplingNn.md)
  : Sampling pseudo-absences for the training and testing datasets.
- [`plotDensityLines()`](https://mdmdma.github.io/USE.MCMC/reference/plotDensityLines.md)
  : plotDensityLines enables the plotting of density function as well as
  a trace of a chain made of points in a dataframe on this density
  surface. In addition the species that was used to generate the density
  model can be supplied to verify the performance of the model. If that
  is the case the supplied dataset will be plotted as points instead of
  as a chain. This function has some issues with the check of devtools
- [`plotInGeographicalSpace()`](https://mdmdma.github.io/USE.MCMC/reference/plotInGeographicalSpace.md)
  : Function that plots the geographical location of points onto a
  raster
- [`rastPCA()`](https://mdmdma.github.io/USE.MCMC/reference/rastPCA.md)
  : Principal Component Analysis for Rasters
- [`thresh.inspect()`](https://mdmdma.github.io/USE.MCMC/reference/thresh.inspect.md)
  : Inspect the effect of the kernel threshold parameter on the
  environmental space partitioning
- [`uniformSampling()`](https://mdmdma.github.io/USE.MCMC/reference/uniformSampling.md)
  : Uniform sampling of the environmental space
