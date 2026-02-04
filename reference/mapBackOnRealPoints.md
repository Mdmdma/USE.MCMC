# Map Back on Real Points searches the closest point in the dataset regarding the given point and the dimensions given

Map Back on Real Points searches the closest point in the dataset
regarding the given point and the dimensions given

## Usage

``` r
mapBackOnRealPoints(dataset = NULL, point = NULL, dim = "", threshold = 0.5)
```

## Arguments

- dataset:

  dataframe of the target dataset

- point:

  dataframe containing the given point

- dim:

  string vector containing the names of the dimensions that should be
  respected

- threshold:

  threshold obove which we consider the distance too big and want to
  discard the point

## Value

closest point that we could find in the dataset with the distance added
as a parameter
