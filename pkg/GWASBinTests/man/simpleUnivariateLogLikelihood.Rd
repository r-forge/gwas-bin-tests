\name{simpleUnivariateLogLikelihood}
\alias{simpleUnivariateLogLikelihood}
\title{compute the univariate log likelihood of a dataset...}
\usage{simpleUnivariateLogLikelihood(genotypes_list, variables, error_models, regCoeff=0)
}
\description{compute the univariate log likelihood of a dataset}
\details{This is a reference implementation, very naive on purpose, optimized for clarity.}
\value{the log-likelihood for the set of markers}
\arguments{\item{genotypes_list}{a list of vector of integers: see \code{\link{fetchGenotypes}}}
\item{variables}{a data frame of variables, possibly empty}
\item{error_models}{list of genotyping error model for each marker, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}
if it is not a list, it must be a matrix and this error model will be used for all markers in the dataset}
\item{regCoeff}{see \code{\link{simpleUnivariateLogLikelihoodForMarker}}}
}

