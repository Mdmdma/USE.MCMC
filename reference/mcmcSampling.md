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
  covariance.correction = 1
)
```

## Arguments

- dataset:

  sf dataframe from which the points are sampled

- dimensions:

  string vector containing the dimensions that should be included in the
  random walk

- densityFunction:

  Function that can take a point given as a sf dataframe as a input and
  returns the target density at that location.

- proposalFunction:

  Function that can take a point given as a sf dataframe and a vector of
  strings specifying the row names that should be changed as a input and
  returns a new proposed point

- n.sample.points:

  Number of points to be sampled

- burnIn:

  Integer, sets the number of samples per adaptive burn in step. If set
  to 0, burn in is skipped

- verbose:

  Boolean to toggle progress updates

- covariance.correction:

  Integer, initial value of the covariance correction.

## Value

A sf dataframe containing the sampled points
