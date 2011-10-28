\name{missingRateFromGenotypes}
\alias{missingRateFromGenotypes}
\title{compute the missing rate from the genotypes of a marker...}
\usage{missingRateFromGenotypes(genotypes)
}
\description{compute the missing rate from the genotypes of a marker}
\details{The missing rate (\eqn{\beta}) is the proportion of missing data (i.e. NA)
among the genotypes}
\value{the missing rate}
\arguments{\item{genotypes}{a vector of integers: see \code{\link{convertGenotypes}}}
}
\examples{data(srdta)
beta <- missingRateFromGenotypes( fetchGenotypes(srdta, 1) )}

