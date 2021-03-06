% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate.msam}
\alias{simulate.msam}
\title{Simulate from a msam object}
\usage{
\method{simulate}{msam}(object, nsim = 1, seed = NULL,
  full.data = TRUE, ready.to.fit = FALSE, ...)
}
\arguments{
\item{object}{msam object result from multisam.fit}

\item{nsim}{Number of simulations}

\item{seed}{random number seed}

\item{full.data}{should a full data set for sam.fit be returned?}

\item{...}{Other arguments not used}
}
\value{
a list of lists.
}
\description{
Simulate from a msam object
}
\author{
Christoffer Moesgaard Albertsen
}
