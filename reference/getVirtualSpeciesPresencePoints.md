# Function to streamline the generation and sampling of a virtual species

Function to streamline the generation and sampling of a virtual species

## Usage

``` r
getVirtualSpeciesPresencePoints(env.data = NULL, n.samples = 0, plot = FALSE)
```

## Arguments

- env.data:

  Raster containing environmental parameters

- n.samples:

  Number of samples

- plot:

  Boolean sets if the function should provide plots of the sampled
  points and the species suitability.

## Value

List containing the sampled points as a dataframe as well as other nice
things that can later be used to plot
