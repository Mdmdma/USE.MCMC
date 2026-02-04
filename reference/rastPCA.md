# Principal Component Analysis for Rasters

The `rastPCA` function calculates the principal component analysis (PCA)
for SpatRaster, RasterBrick, or RasterStack objects and returns a
SpatRaster with multiple layers representing the PCA components.
Internally, `rastPCA` utilizes the
[princomp](https://rdrr.io/r/stats/princomp.html) function for R-mode
PCA analysis. The covariance matrix is computed using all the
observations within the provided SpatRaster object, which describes the
environmental conditions. The covariance matrix obtained is subsequently
utilized as input for the `princomp` function, which conducts the PCA.
The resulting PCA components are then used to generate the final
SpatRaster, consisting of multiple layers that represent the PCA
components.

## Usage

``` r
rastPCA(env.rast, nPC = NULL, naMask = TRUE, stand = FALSE)
```

## Arguments

- env.rast:

  A RasterStack, RasterBrick or a SpatRaster object comprising the
  variables describing the environmental space.

- nPC:

  Integer. Number of PCA components to return.

- naMask:

  Logical. Masks all pixels which have at least one NA (default `TRUE`
  is recommended but introduces a slow-down.

- stand:

  Logical. If `TRUE`, perform standardized PCA. Corresponds to centered
  and scaled input image. This is usually beneficial for equal weighting
  of all layers. (`FALSE` by default)

## Value

Returns a named list containing the PCA model object (\$pca) and the
SpatRaster with the principal component layers (\$PCs).

## Details

Pixels with missing values in one or more bands will be set to NA. The
built-in check for such pixels can lead to a slow-down of rastPCA.
However, if you make sure or know beforehand that all pixels have either
only valid values or only NAs throughout all layers you can disable this
check by setting `naMask=FALSE` which speeds up the computation.

Standardized PCA (`stand=TRUE`) can be useful if imagery or bands of
different dynamic ranges are combined. In this case, the correlation
matrix is computed instead of the covariance matrix, which has the same
effect as using normalised bands of unit variance.

## See also

The `rastPCA` function has been conceptualized starting from
`RStoolbox::rasterPCA` (<https://github.com/bleutner/RStoolbox>).
