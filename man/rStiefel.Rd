% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{rStiefel}
\alias{rStiefel}
\title{Random sample from Stiefel manifold.}
\usage{
rStiefel(p, q)
}
\arguments{
\item{p}{row dimension}

\item{q}{col dimension}
}
\value{
A \eqn{p \times q}{p x q} semi-orthogonal matrix.
}
\description{
Draws a random sample from the invariant measure on the Stiefel manifold
\eqn{S(p, q)}.
}
\examples{
V <- rStiefel(6, 4)
}
