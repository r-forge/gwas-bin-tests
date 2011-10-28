\name{simpleDivariatePvalue}
\alias{simpleDivariatePvalue}
\title{compute the divariate likelihood Ratio pvalue for a dataset/bin using permutations...}
\usage{simpleDivariatePvalue(nb_permutations, genotypes_list, phenotypes, error_models,
    covariables=data.frame(), regCoeff=0)
}
\description{compute the divariate likelihood Ratio pvalue for a dataset/bin using permutations}
\details{See \code{\link{simpleDivariateLogLikelihood}}
This is a reference implementation, very naive on purpose, optimized for clarity.

If there is only one SNP or less, it stops with an error message!}
\value{a list (pvalue=, observed_score=, nb_permutations=)}
\arguments{\item{nb_permutations}{the number of permutations to do}
\item{genotypes_list}{a list of vector of integers: see \code{\link{fetchGenotypes}}}
\item{phenotypes}{a vector of phenotypes of same length as the elements of genotypes_list}
\item{error_models}{list of genotyping error model for each marker, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}
if it is not a list, it must be a matrix and this error model will be used for all markers in the dataset}
\item{covariables}{a data frame of cvariables: \eqn{ (\mathcal{V})_m }}
\item{regCoeff}{see \code{\link{simpleDivariateLogLikelihoodForPair}}}
}

