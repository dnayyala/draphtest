% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dir_hessian.R
\name{dir.hessian}
\alias{dir.hessian}
\title{Hessian for Dirichlet distribution}
\usage{
dir.hessian(x, mu, alpha, param = "alpha")
}
\arguments{
\item{x}{The data matrix of size \eqn{p \times n} which is generated from \eqn{p}-dimensional Dirichlet distribution.}

\item{mu}{A positive scalar parameter which represents the dispersion parameter of the Dirichlet distribution.}

\item{alpha}{A \eqn{p \times 1} parameter vector which is the log-transformed mean of Dirichlet distribution.}

\item{param}{A character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".}
}
\value{
the Hessian matrix of Dirichlet distribution size of \eqn{(p-1) \times (p-1)}
Off-diagonal elements of Hessian Matrix
Diagonal elements of Hessian Matrix
}
\description{
Hessian matrix of Dirichlet distribution with respect to the parameters of the distribution: alpha and mu
}
