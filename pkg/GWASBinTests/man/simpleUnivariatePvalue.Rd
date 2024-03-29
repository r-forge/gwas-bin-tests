\name{simpleUnivariatePvalue}
\alias{simpleUnivariatePvalue}
\title{compute the univariate likelihood Ratio pvalue for a dataset/bin using permutations...}
\usage{simpleUnivariatePvalue(nb_permutations, genotypes_list, phenotypes, error_models,
    covariables=data.frame(), regCoeff=0)
}
\description{compute the univariate likelihood Ratio pvalue for a dataset/bin using permutations}
\details{See \code{\link{simpleUnivariateLogLikelihoodRatio}}
This is a reference implementation, very naive on purpose, optimized for clarity.}
\value{a list (pvalue=, observed_score=, nb_permutations=)}
\arguments{\item{nb_permutations}{the number of permutations to do}
\item{genotypes_list}{a list of vector of integers: see \code{\link{fetchGenotypes}}}
\item{phenotypes}{a vector of phenotypes of same length as the elements of genotypes_list}
\item{error_models}{list of genotyping error model for each marker, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}
if it is not a list, it must be a matrix and this error model will be used for all markers in the dataset}
\item{covariables}{a data frame of cvariables: \eqn{ (\mathcal{V})_m }}
\item{regCoeff}{see \code{\link{simpleUnivariateLogLikelihoodForMarker}}}
}

