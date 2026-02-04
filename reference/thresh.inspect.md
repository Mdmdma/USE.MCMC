# Inspect the effect of the kernel threshold parameter on the environmental space partitioning

`thresh.inspect` function allows for a pre-inspection of the impact that
selecting a specific threshold for the kernel-based filter will have on
the exclusion of the environmental space in the subsequent uniform
sampling of the pseudo-absences process (see `paSampling`). By providing
a range of threshold values, the function generates a plot that
illustrates the entire environmental space, including the portion
delineated by the kernel-based filter and the associated convex-hull.
This plot helps visualize the areas that will be excluded from the
uniform sampling of the pseudo-absences. This functionality proves
particularly valuable in determining a meaningful threshold for the
kernel-based filter in specific ecological scenarios. For instance, when
dealing with sink populations, selecting the appropriate threshold
enables the exclusion of environmental space regions where the species
is present, but the conditions are unsuitable. This allows for a more
accurate sampling of pseudo-absences, considering the unique
requirements of different ecological contexts.

## Usage

``` r
thresh.inspect(env.rast, pres = NULL, thres = 0.75, H = NULL)
```

## Arguments

- env.rast:

  A RasterStack, RasterBrick or a SpatRaster object comprising the
  variables describing the environmental space.

- pres:

  A SpatialPointsDataframe, a SpatVector or an sf object including the
  presence-only observations of the species of interest.

- thres:

  (double) This value or vector of values identifies the quantile value
  used to specify the boundary of the kernel density estimate (default
  `thres=0.75` ). Thus, probability values higher than the threshold
  should indicate portions of the multivariate space likely associated
  with presence points.

- H:

  The kernel bandwidth (i.e., the width of the kernel density function
  that defines its shape) excluding the portion of the environmental
  space associated with environmental conditions likely suitable for the
  species. It can be either defined by the user or automatically
  estimated by `paSampling` via
  [`ks::Hpi`](https://rdrr.io/pkg/ks/man/Hpi.html).

## Value

A ggplot2 object showing how the environmental space is partitioned
accordingly to the selected `thres` values.
