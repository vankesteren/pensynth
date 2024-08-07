% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{rarnorm}
\alias{rarnorm}
\title{Generate data from normal distribution with AR1 parameter}
\usage{
rarnorm(n, mean = 0, sd = 1, phi = 0)
}
\arguments{
\item{n}{number of observations.}

\item{mean}{marginal mean}

\item{sd}{marginal standard deviation}

\item{phi}{autoregressive parameter (-1 < phi < 1)}
}
\value{
vector of numeric values
}
\description{
Generate data from normal distribution with AR1 parameter
}
\details{
Note that, unlike \code{\link[stats:Normal]{stats::rnorm()}}, this function is
not vectorized over mean, sd, or phi.
}
\keyword{internal}
