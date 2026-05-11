# Batched GMM density evaluation using precomputed parameters

Vectorized variant of `fast_gmm_density` for evaluating many points at
once. Operates per component on the whole batch of points, then combines
via log-sum-exp. Substantially faster than looping single-point
evaluations when the caller already has all query points in hand (e.g.
plotting density rasters, threshold sweeps, or sweep-style diagnostics
in the vignettes).

## Usage

``` r
fast_gmm_density_batch(X, pre)
```

## Arguments

- X:

  Numeric matrix (n x d) or data.frame coercible to one — query points
  in rows, dimensions in columns.

- pre:

  List of precomputed parameters from `precompute_gmm_params`.

## Value

Numeric vector of length n with the density at each row of X.
