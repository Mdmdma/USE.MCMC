# MCMC sampling from a given dataset

MCMC sampling from a given dataset

## Usage

``` r
mcmcSampling(
  dataset = NULL,
  dimensions = list(""),
  densityFunction = alwaysOne,
  proposalFunction = addHighDimGaussian(dim = length(dimensions)),
  n.sample.points = 0,
  burnIn = 1000,
  verbose = TRUE,
  covariance.correction = 1,
  max.burnin.cycles = 50,
  engine = c("auto", "R", "cpp")
)
```

## Arguments

- dataset:

  sf dataframe from which the points are sampled

- dimensions:

  string vector containing the dimensions that should be included in the
  random walk

- densityFunction:

  Function that can take a point given as a numeric vector as input and
  returns the target density at that location.

- proposalFunction:

  Function that can take a point given as a numeric vector and a
  covariance adjuster as input and returns a new proposed point as a
  numeric vector.

- n.sample.points:

  Number of points to be sampled

- burnIn:

  Integer, number of Robbins-Monro burn-in adaptation steps performed
  before sampling. During each step the proposal scale is adjusted
  toward target acceptance 0.234 (Roberts/Rosenthal 2009). Set to 0 to
  skip adaptation and start sampling immediately at the user-supplied
  `covariance.correction`.

- verbose:

  Boolean to toggle progress updates

- covariance.correction:

  Integer, initial value of the covariance correction.

- max.burnin.cycles:

  Deprecated. Retained for backwards compatibility; ignored by the
  current Robbins-Monro burn-in.

- engine:

  One of `"auto"` (default), `"R"`, or `"cpp"`. `"auto"` picks the C++
  inner loop when both `densityFunction` and `proposalFunction` are
  built by
  [`mclustDensityFunction()`](https://mdmdma.github.io/USE.MCMC/reference/mclustDensityFunction.md)
  and
  [`addHighDimGaussian()`](https://mdmdma.github.io/USE.MCMC/reference/addHighDimGaussian.md)
  (they carry the required `rcpp_spec` attribute) and falls back to the
  R loop otherwise. `"cpp"` forces the C++ path and errors if a custom
  closure is supplied. `"R"` forces the pure-R reference loop.

## Value

A data.frame containing the sampled points with dimension columns and a
density column
