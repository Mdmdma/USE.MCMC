# Function that plots the geographical location of points onto a raster

Function that plots the geographical location of points onto a raster

## Usage

``` r
plotInGeographicalSpace(
  presence.distribution.raster = NULL,
  presence.points = NULL,
  absence.points = NULL,
  minimal = FALSE
)
```

## Arguments

- presence.distribution.raster:

  Raster containing the target distribution

- presence.points:

  Dataframe or list containing a geometry cloumn that contains points to
  be plotted

- absence.points:

  Dataframe or list containing a geometry cloumn that contains points to
  be plotted

- minimal:

  Boolean if TRUE removes the labels legend and title

## Value

NULL, as the function just plots
