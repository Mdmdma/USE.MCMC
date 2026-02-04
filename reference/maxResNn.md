# maxResNN is a function that can be used to compute a reasonable grid resolution for nearest neighbor based uniform sampling. The core idea behind its working principle is that we want to expect a grid cell to contain points if it overlaps with the environment. This implementation looks at the low density regions using distance to n neighbors as a proxy. The approach assumes that both coordinates have similar range, as the axes are not weighted when computing Neighbors and converting from distances to number of grid cells If multiple points from the same grid cell should be sampled, the number of neighbors included in the computation should be set accordingly

maxResNN is a function that can be used to compute a reasonable grid
resolution for nearest neighbor based uniform sampling. The core idea
behind its working principle is that we want to expect a grid cell to
contain points if it overlaps with the environment. This implementation
looks at the low density regions using distance to n neighbors as a
proxy. The approach assumes that both coordinates have similar range, as
the axes are not weighted when computing Neighbors and converting from
distances to number of grid cells If multiple points from the same grid
cell should be sampled, the number of neighbors included in the
computation should be set accordingly

## Usage

``` r
maxResNn(
  env.data.raster,
  dimensions = c("PC1", "PC2"),
  low.end.of.inclueded.points = 20,
  high.end.of.included.points = 4,
  n.neighbors = 2,
  PCA = FALSE
)
```

## Arguments

- env.data.raster:

  Raster containing the environmental parameters

- dimensions:

  vector containing the dimensions that should be used for the grid
  computation

- low.end.of.inclueded.points:

  low cutoff

- high.end.of.included.points:

  high cutoff

- n.neighbors:

  number of neighbors used for the computation

- PCA:

  can be set to true if rastPCA was already used to perform a pca to
  save time recomputing

## Value

maximal number of grid cells useful for the nearest neighbor based
approach
