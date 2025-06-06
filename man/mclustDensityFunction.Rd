% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mclustDensityFunction.R
\name{mclustDensityFunction}
\alias{mclustDensityFunction}
\title{Helper to create a Density function that uses mclust Gaussian mixtures}
\usage{
mclustDensityFunction(
  env.model = NULL,
  species.model = NULL,
  dim = "",
  threshold = 0.01,
  species.cutoff.threshold = 0.1
)
}
\arguments{
\item{env.model}{mclust gaussian mixture that uses points}

\item{species.model}{mclust gaussian mixture that uses points}

\item{dim}{string vector specifing the names of the dimensions}

\item{threshold}{sets the curoff density from the environment}

\item{species.cutoff.threshold}{set the scaling factor with which the species model gets scaled before subtraction.}
}
\value{
Function that can calculates the density at a point
}
\description{
Helper to create a Density function that uses mclust Gaussian mixtures
}
