\name{convertGenotypes}
\alias{convertGenotypes}
\title{convert and check the genotypes to our simple pure-R implementation conventions...}
\usage{convertGenotypes(genotypes)
}
\description{convert and check the genotypes to our simple pure-R implementation conventions}
\details{The genotypes must be an integer vector.
The code for missing data is -1, so all NAs will be converted into -1

All other values will raise an error.

This format of genotypes is expected by all the simple_* functions.
To extract genotypes from a Gws, rather use the \code{\link{fetchGenotypes}} and \code{\link{fetchGenotypesAsList}}
functions that will call this function anyway.}
\value{the cleaned and checked integer vector of genotypes, with values in -1:2}
\arguments{\item{genotypes}{a vector (or matrix) of integers
with values meaning: \describe{
\item{homozygous1 (AA) }{0}
\item{heterozygous (AB) }{1}
\item{homozygous2 (BB) }{2}
\item{missing (\eqn{\emptyset}) }{\code{NA} or -1}
}. There is no assumption about which allele is minor.}
}
\examples{data(srdta)
gs <- convertGenotypes( as.integer( as.double(srdta[,1:20]) ) )}

