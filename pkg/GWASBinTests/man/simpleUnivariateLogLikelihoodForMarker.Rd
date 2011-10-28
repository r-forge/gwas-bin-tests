\name{simpleUnivariateLogLikelihoodForMarker}
\alias{simpleUnivariateLogLikelihoodForMarker}
\title{compute the univariate log likelihood for a given marker...}
\usage{simpleUnivariateLogLikelihoodForMarker(genotypes, variables, error_model, regCoeff=0)
}
\description{compute the univariate log likelihood for a given marker}
\details{This is the a very simple implementation, self-contained and easy to check
but not performant. It is used as illustration or to check other implementations.}
\value{the log likelihood numeric value. If there is no genotype for a category, the likelihood will not be defined
and so it will return NaN}
\arguments{\item{genotypes}{a vector of integers: see \code{\link{convertGenotypes}}}
\item{variables}{a data frame of variables, possibly empty}
\item{error_model}{a genotyping error model, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}}
\item{regCoeff}{regularization rate: will be used to compute the constant to add
in each cell of the contingency table to regularize estimators. This constant will
be = rate*total/nb_of_cells. See \strong{\link{concept#regularization.of.table.estimates}}.}
}

