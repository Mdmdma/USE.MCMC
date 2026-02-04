# Virtual species probability of occurrence

The `SpatialProba` function calculates the simulated probability of
occurrence of a virtual species based on an additive model that
incorporates environmental variables. The model considers both linear
and quadratic relationships between the environmental factors and the
species' probability of presence. This function uses environmental data
provided as a SpatRaster object (e.g., temperature, precipitation) to
compute the probability of species presence across a defined area of
interest. The resulting probabilities are mapped to a range between 0
and 1, representing the likelihood of species occurrence in the given
locations.

## Usage

``` r
SpatialProba(coefs, env.rast, quadr_term, marginalPlots)
```

## Arguments

- coefs:

  a named vector of regression parameters. Names must match those of the
  environmental layers (except for intercept, and quadratic terms).
  Parameters for quadratic terms must have the prefix 'quadr\_' (e.g.,
  `quadr_bio1`).

- env.rast:

  a SpatRaster object with environmental layers to generate the spatial
  layer of probabilities.

- quadr_term:

  a named vector with names of coefs for which a quadratic term is
  specified (without prefix 'quadr\_').

- marginalPlots:

  logical, if TRUE, returns marginal plots.
