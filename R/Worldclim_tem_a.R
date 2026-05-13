#' WorldClim bioclimatic variables from the original USE package
#'
#' A subset of WorldClim bioclimatic variables cropped on Central and Western
#' Europe. Identical to \code{USE::Worldclim_tmp} from the upstream
#' \href{https://github.com/danddr/USE}{USE} package, included here so this
#' package's examples and tests can reproduce the upstream demos without
#' requiring \pkg{USE} to be installed.
#'
#' @docType data
#' @keywords datasets
#' @name Worldclim_tem_a
#' @usage data(Worldclim_tem_a)
#' @format A data frame with 18953 rows and 8 columns: \code{x}, \code{y},
#'   and the bioclim layers \code{bio1}, \code{bio3}, \code{bio9},
#'   \code{bio12}, \code{bio13}, \code{bio15}.
#' @source \code{geodata::worldclim_global(var='bio', res=10, path=getwd())[[c(1, 3, 9, 12, 13, 15)]]}
'Worldclim_tem_a'
