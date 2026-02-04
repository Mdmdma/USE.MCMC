# Helper to create a Density function that uses mclust Gaussian mixtures

Helper to create a Density function that uses mclust Gaussian mixtures

## Usage

``` r
mclustDensityFunction(
  env.model = NULL,
  species.model = NULL,
  dim = "",
  threshold = 0.01,
  species.cutoff.threshold = 0.1
)
```

## Arguments

- env.model:

  mclust gaussian mixture that uses points

- species.model:

  mclust gaussian mixture that uses points

- dim:

  string vector specifing the names of the dimensions

- threshold:

  sets the curoff density from the environment

- species.cutoff.threshold:

  set the scaling factor with which the species model gets scaled before
  subtraction.

## Value

Function that can calculates the density at a point
