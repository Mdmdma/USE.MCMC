% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highDimGaussian.R
\name{highDimGaussian}
\alias{highDimGaussian}
\title{multivariate gaussian random number generator}
\usage{
highDimGaussian(dim = 1, mean_vec = matrix(1, ncol = dim), cov_mat = diag(dim))
}
\arguments{
\item{dim}{number of dimensions}

\item{mean_vec}{first moment}

\item{cov_mat}{second moment}
}
\value{
a function that generates samples from a higher dimensional gaussian
}
\description{
multivariate gaussian random number generator
}
\keyword{internal}
