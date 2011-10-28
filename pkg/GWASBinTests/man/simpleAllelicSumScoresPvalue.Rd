\name{simpleAllelicSumScoresPvalue}
\alias{simpleAllelicSumScoresPvalue}
\title{compute the sum of allelic pearson score pvalue on a dataset using permutations...}
\usage{simpleAllelicSumScoresPvalue(nb_permutations, genotypes_list, phenotypes, covariables=data.frame(),
    regCoeff=0)
}
\description{compute the sum of allelic pearson score pvalue on a dataset using permutations}
\details{See \code{\link{simpleAllelicSumScores}}
This is a reference implementation, very naive on purpose, optimized for clarity.}
\value{a list (pvalue=, observed_score=, nb_permutations=)}
\arguments{\item{nb_permutations}{the number of permutations to do}
\item{genotypes_list}{a list of vector of integers: see \code{\link{fetchGenotypesAsList}}}
\item{phenotypes}{a vector of phenotypes of same length as the elements of genotypes_list}
\item{covariables}{a data frame of cvariables: \eqn{ (\mathcal{V})_m }}
\item{regCoeff}{regularization rate: will be used to compute the constant to add
in each cell of the contingency table to regularize estimators. This constant will
be = rate*total/nb_of_cells. See \strong{\link{concept#regularization.of.table.estimates}}.}
}

