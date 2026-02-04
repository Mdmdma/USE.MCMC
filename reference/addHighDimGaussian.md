# addHighDimGaussian

Adds high dimensional Gaussian noise to specific parameters of a point
given as a dataframe.

## Usage

``` r
addHighDimGaussian(
  dim = 0,
  mean.vec = matrix(0, ncol = dim),
  cov.mat = diag(dim)
)
```

## Arguments

- dim:

  integer that specifies the number of dimensions that should be altered

- mean.vec:

  vector of the means of the Gaussian to be added

- cov.mat:

  covariance matrix of the Gaussian to be added

## Value

function that takes a point given as a dataframe as input and returns it
with Gaussian noise added to the specified dimensions
