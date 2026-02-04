# Uniform sampling of the environmental space

`uniformSampling` performs the uniform sampling of observations within
the environmental space. Note that `uniformSampling` can be more
generally used to sample observations (not necessarily associated with
species occurrence data) within bi-dimensional spaces (e.g., vegetation
plots). Being designed with species distribution models in mind,
`uniformSampling` allows collectively sampling observations for both the
training and testing dataset (optional). In both cases, the user must
provide a number of observations that will be sampled in each cell of
the sampling grid (`n.tr`: points for the training dataset; `n.ts`:
points for the testing dataset). Note that the optimal resolution of the
sampling grid can be found using the `optimRes` function.

## Usage

``` r
uniformSampling(
  sdf,
  grid.res,
  n.tr = 5,
  n.prev = NULL,
  sub.ts = FALSE,
  n.ts = 5,
  plot_proc = FALSE,
  verbose = FALSE
)
```

## Arguments

- sdf:

  an sf object having point geometry given by the PC-scores values

- grid.res:

  (integer) resolution of the sampling grid. The resolution can be
  arbitrarily selected or defined using the
  [`optimRes()`](https://mdmdma.github.io/USE.MCMC/reference/optimRes.md)
  function.

- n.tr:

  (integer; optional) number of expected points given a certain
  prevalence threshold for the training dataset.

- n.prev:

  (double) sample prevalence

- sub.ts:

  (logical) sample the validation points

- n.ts:

  (integer; optional) number of points for the testing dataset to sample
  in each cell of the sampling grid. sub.ts argument must be TRUE.

- plot_proc:

  (logical) plot progress of the sampling

- verbose:

  (logical) Print verbose

## Value

An sf object with the coordinates of the sampled points both in the
geographical and environmental space
