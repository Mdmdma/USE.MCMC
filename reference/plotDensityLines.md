# plotDensityLines enables the plotting of density function as well as a trace of a chain made of points in a dataframe on this density surface. In addition the species that was used to generate the density model can be supplied to verify the performance of the model. If that is the case the supplied dataset will be plotted as points instead of as a chain. This function has some issues with the check of devtools

plotDensityLines enables the plotting of density function as well as a
trace of a chain made of points in a dataframe on this density surface.
In addition the species that was used to generate the density model can
be supplied to verify the performance of the model. If that is the case
the supplied dataset will be plotted as points instead of as a chain.
This function has some issues with the check of devtools

## Usage

``` r
plotDensityLines(
  dataset,
  xlim = c(0, 1),
  ylim = c(0, 1),
  title = "Connected Data Points",
  cols = NULL,
  lines = FALSE,
  density = FALSE,
  species = NULL,
  densityFunction = alwaysOne,
  resolution = 10,
  minimal = FALSE
)
```

## Arguments

- dataset:

  Dataframe specifying the dataset to be plotted as a line or points
  respectively

- xlim:

  x-limit of the plot

- ylim:

  y-limit of the plot

- title:

  title of the plot

- cols:

  specifies the columns of the dataframe that are used for the density
  computation. It has to match to the columns used to build the density
  model.

- lines:

  boolean that tells the function if the lines should be plotted. If
  false the dataset will be ploted as points

- density:

  boolean that tells the function if the density should be plotted

- species:

  dataframe containing the simulated species, it has to contain the
  column names of cols

- densityFunction:

  a function that can take a dataframe with the columns given in cols
  and retruns the density at that point

- resolution:

  int sets the resolution of the denisity grid

- minimal:

  boolean if true removes titel, label, legend and axis

## Value

the greated plot
