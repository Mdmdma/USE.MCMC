# addHighDimGaussian

Adds high dimensional Gaussian noise to specific parameters of a point
given as a numeric vector.

## Usage

``` r
addHighDimGaussian(
  dim = 0,
  mean.vec = matrix(0, ncol = dim),
  cov.mat = diag(dim),
  batch.size = 1000
)
```

## Arguments

- dim:

  integer that specifies the number of dimensions that should be altered

- mean.vec:

  vector of the means of the Gaussian to be added

- cov.mat:

  covariance matrix of the Gaussian to be added

- batch.size:

  integer, number of random vectors to pre-generate at a time for
  efficiency

## Value

function that takes a point given as a numeric vector and returns it
with Gaussian noise added
