% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/acceptNextPoint.R
\name{acceptNextPoint}
\alias{acceptNextPoint}
\title{Accept next point ?}
\usage{
acceptNextPoint(current.point, proposed.point)
}
\arguments{
\item{current.point}{Current state of the markov chain, given as a sf dataframe with a column called $density containing the value of the target density.}

\item{proposed.point}{Proposed next state of the markov chain, given as a sf dataframe with a column called $density containing the value of the target density.}
}
\value{
TRUE if the proposed point should be accepted
}
\description{
\code{acceptNextPoint} checks whether the proposed point should be accepted
}
\keyword{internal}
