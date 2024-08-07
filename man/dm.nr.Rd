% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dm_nr.R
\name{dm.nr}
\alias{dm.nr}
\title{Newton-Raphson for Dirichlet-multinomial}
\usage{
dm.nr(x, mu = NULL, alpha = NULL, tol = 1e-04, maxiters = 100)
}
\arguments{
\item{x}{A data matrix of size \eqn{p \times n} generated from \link{p}-dimensional Dirichlet-multinomial distribution.}

\item{mu}{A scalar parameter which represents the dispersion parameter of the distribution}

\item{alpha}{Vector of length p which represents the log-transformed mean vector of the distribution.}

\item{tol}{Tolerance cutoff for determining convergence of the Newton-Raphson algorithm}

\item{maxiters}{Maximum number of iterations for the algorithm. Default is set at 100.}
}
\value{
the function will return a vector of length \eqn{p} containing the estimate of the param alpha and a the estimate of mu
}
\description{
Newton-Rapshon algorithm to estimate parameter alpha (or theta) of the Dirichlet-multinomial distribution
}
