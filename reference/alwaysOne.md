# AlwaysOne

As the name implies, this function returns always 1, no matter the
input. It can be used as a dummy desity function to perform uniform
sampling of higher dimensional space. As it works with any Input, it can
also be used as a placeholder to evade errors. always

## Usage

``` r
alwaysOne(...)
```

## Arguments

- ...:

  Any input is possible, as it does not have any effect

## Value

1, always

## Examples

``` r
alwaysOne(stop("even this input is fine"))
#> [1] 1

```
