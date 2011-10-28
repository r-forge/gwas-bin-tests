\name{simpleDivariateLogLikelihoodForPair}
\alias{simpleDivariateLogLikelihoodForPair}
\title{compute the divariate log likelihood for a given pair of markers...}
\usage{simpleDivariateLogLikelihoodForPair(g1, e1, g2, e2, variables, regCoeff=0)
}
\description{compute the divariate log likelihood for a given pair of markers}
\details{This is the a very simple implementation, self-contained and easy to check
but not performant. It is used as illustration or to check other implementations.}
\value{the log likelihood numeric value. If there is no genotype for a category, the likelihood will not be defined
and so it will return NaN}
\arguments{\item{g1}{the genotypes of the first marker as a vector of integers: see \code{\link{convertGenotypes}}}
\item{e1}{the genotyping error model for the first marker, see \link{concept#genotype.error.model},
\code{\link{simpleGenotypeErrorModel}}}
\item{g2}{the genotypes of the second marker}
\item{e2}{the genotyping error model for the second marker}
\item{variables}{a data frame of variables, possibly empty}
\item{regCoeff}{regularization rate: will be used to compute the constant to add
in each cell of the contingency table to regularize estimators. This constant will
be = rate*total/nb_of_cells. See \strong{\link{concept#regularization.of.table.estimates}}.}
}

