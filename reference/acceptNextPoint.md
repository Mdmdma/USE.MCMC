# Accept next point ?

`acceptNextPoint` checks whether the proposed point should be accepted

## Usage

``` r
acceptNextPoint(current.point, proposed.point)
```

## Arguments

- current.point:

  Current state of the markov chain, given as a sf dataframe with a
  column called \$density containing the value of the target density.

- proposed.point:

  Proposed next state of the markov chain, given as a sf dataframe with
  a column called \$density containing the value of the target density.

## Value

TRUE if the proposed point should be accepted
