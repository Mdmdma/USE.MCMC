#!/usr/bin/env bash
#SBATCH --job-name=use-mcmc-nn
#SBATCH --time=00:50:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/scratch/%u/use-mcmc-nn-%j.out

# Run the USE.MCMC testthat suite + regenerate roxygen docs on a compute node, inside
# the rocker-rstudio R 4.5 apptainer container. The container/base image lacks the
# GDAL/PROJ/GEOS/udunits system libs that terra/sf link against, so we inject a
# soname-matching set via LD_LIBRARY_PATH: GDAL+PROJ (libgdal.so.33 / libproj.so.25)
# from the manual paraview 6.0.1 install, GEOS+udunits from the 2025-06 module stack.
# (The tests do not read/write geospatial files, so GDAL drivers are not exercised --
# terra/sf only need to load.) Compiling the package src must not run on the login node.

set -euo pipefail

SIF=/cluster/scratch/${USER}/rocker_rstudio_4.5.sif
PKG=/cluster/home/${USER}/USE.MCMC
RLIB=/cluster/home/${USER}/R/rocker-rstudio/4.5

PV=/cluster/software/manual/paraview/6.0.1/x86_64/lib
GEOS=/cluster/software/stacks/2025-06/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0/geos-3.13.0-5e4em54xy2jj543vakxivg7fnblmjm6b/lib
UDU=/cluster/software/stacks/2025-06/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0/udunits-2.2.28-xas2wcn64didtblh7nbfypjox5b2rdqm/lib

cd "$PKG"

export APPTAINERENV_R_LIBS="$RLIB"
export APPTAINERENV_LD_LIBRARY_PATH="$PV:$GEOS:$UDU:/usr/local/lib/R/lib:/usr/lib/x86_64-linux-gnu:/.singularity.d/libs"

apptainer exec --bind /cluster "$SIF" Rscript -e '
  .libPaths(c("'"$RLIB"'", .libPaths()))
  options(crayon.enabled = FALSE, warn = 1)
  cat("R:", R.version.string, "\n")
  cat("terra loads:", requireNamespace("terra", quietly = TRUE),
      "| sf loads:", requireNamespace("sf", quietly = TRUE), "\n")

  cat("\n== devtools::test() ==\n")
  tests_ok <- TRUE
  res <- tryCatch(
    as.data.frame(devtools::test(reporter = "summary", stop_on_failure = FALSE)),
    error = function(e) { cat("TEST RUN ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) {
    tests_ok <- FALSE
  } else {
    bad <- res[res$failed > 0 | res$error, , drop = FALSE]
    cat(sprintf("\nPASS=%d  FAIL=%d  WARN=%d  SKIP=%d\n",
                sum(res$passed), sum(res$failed), sum(res$warning), sum(res$skipped)))
    # focus report on the three NN files we changed
    nn <- res[grepl("Nn", res$file), c("file","test","passed","failed","warning")]
    cat("\n-- NN test files --\n"); print(nn)
    if (nrow(bad) > 0) {
      cat("\n!! FAILURES / ERRORS:\n"); print(bad[, c("file","test","failed","error")]); tests_ok <- FALSE
    }
  }

  cat("\n== roxygen2::roxygenise() ==\n")
  docs_ok <- tryCatch({ roxygen2::roxygenise(); TRUE },
                      error = function(e) { cat("DOC ERROR:", conditionMessage(e), "\n"); FALSE })

  cat(sprintf("\n== RESULT: tests_ok=%s  docs_ok=%s ==\n", tests_ok, docs_ok))
  if (!isTRUE(tests_ok)) quit(status = 1)
'
