# Precompute GMM parameters for fast density evaluation

Extracts and precomputes inverse covariance matrices and log
normalization constants from a densityMclust model to avoid recomputing
them on every call.

## Usage

``` r
precompute_gmm_params(model)
```

## Arguments

- model:

  A densityMclust object

## Value

A list with precomputed parameters
