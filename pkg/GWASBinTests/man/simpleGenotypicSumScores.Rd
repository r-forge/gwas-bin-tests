\name{simpleGenotypicSumScores}
\alias{simpleGenotypicSumScores}
\title{compute the sum of genotypic pearson score on a dataset...}
\usage{simpleGenotypicSumScores(genotypes_list, phenotypes, covariables=data.frame(), regCoeff=0)
}
\description{compute the sum of genotypic pearson score on a dataset}
\details{cf \link{concept#genotypic.sum.score}
This is a reference implementation, very naive on purpose, optimized for clarity.}
\value{the score for the set of markers}
\arguments{\item{genotypes_list}{a list of vector of integers: see \code{\link{fetchGenotypesAsList}}}
\item{phenotypes}{a vector of phenotypes of same length as the elements of genotypes_list}
\item{covariables}{a data frame of cvariables: \eqn{ (\mathcal{V})_m }}
\item{regCoeff}{regularization rate: will be used to compute the constant to add
in each cell of the contingency table to regularize estimators. This constant will
be = rate*total/nb_of_cells. See \strong{\link{concept#regularization.of.table.estimates}}.}
}

