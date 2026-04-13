# Fast GMM density evaluation using precomputed parameters

Computes the density of a Gaussian mixture model at a point using
precomputed inverse covariance matrices and log normalization constants.
Uses log-sum-exp for numerical stability.

## Usage

``` r
fast_gmm_density(x, pre)
```

## Arguments

- x:

  Numeric vector, the point at which to evaluate the density.

- pre:

  List of precomputed parameters from `precompute_gmm_params`.

## Value

Numeric scalar, the density at x.
