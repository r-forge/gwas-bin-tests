

#' compute the sum of genotypic pearson score on a dataset
#'
#' cf \link{concept#genotypic.sum.score}
#' This is a reference implementation, very naive on purpose, optimized for clarity.
#'
#' @param genotypes_list a list of vector of integers: see \code{\link{fetchGenotypesAsList}}
#' @param phenotypes a vector of phenotypes of same length as the elements of genotypes_list
#' @param covariables a data frame of cvariables: \eqn{ (\mathcal{V})_m }
#' @param  regCoeff 	regularization rate: will be used to compute the constant to add
#' 	in each cell of the contingency table to regularize estimators. This constant will
#' 	be = rate*total/nb_of_cells. See \strong{\link{concept#regularization.of.table.estimates}}.
#' @return the score for the set of markers
#' @export
simpleGenotypicSumScores <- function(genotypes_list,  phenotypes, covariables=data.frame(), regCoeff=0) {
	nb_markers <- length(genotypes_list)

	if ( ! is.data.frame(covariables) )
		stop("covariables must be a data frame")

	if ( nb_markers < 1)
		stop("there must be at least genotypes for one marker")

	sum( sapply(genotypes_list, function(g) {
		scores <- by(data.frame(g,phenotypes), covariables, function(x) {

			##  === compute the contingency table  ===
			ct <- table(factor(x[,1], levels=c(0:2), labels=c("aa", "Aa", "AA")), factor(x[,2]))

			# only missing data ==> we return a score of zero so that it has no impact on the global sum
			if ( sum(ct) == 0 )
				return(0)

			# regularization constant if any
			if ( regCoeff != 0) {
				const <- regCoeff * sum(ct) / 3 # compute the constant
				ct <- ct + const # add it to each cell
			}

			expected <- ( rowSums(ct) %o% colSums(ct) ) / sum(ct)
			m <- (ct-expected)^2 / expected
			score <- sum(  m[expected != 0 ] )

			return(score)
		})
	}))
}

#' compute the sum of genotypic pearson score pvalue on a dataset using permutations
#'
#' See \link{concept#genotypic.sum.score} and \code{\link{simpleGenotypicSumScores}}
#' This is a reference implementation, very naive on purpose, optimized for clarity.
#'
#' @param nb_permutations the number of permutations to do
#' @param genotypes_list a list of vector of integers: see \code{\link{fetchGenotypesAsList}}
#' @param phenotypes a vector of phenotypes of same length as the elements of genotypes_list
#' @param covariables a data frame of cvariables: \eqn{ (\mathcal{V})_m }
#' @param  regCoeff 	regularization rate: will be used to compute the constant to add
#' 	in each cell of the contingency table to regularize estimators. This constant will
#' 	be = rate*total/nb_of_cells. See \strong{\link{concept#regularization.of.table.estimates}}.
#'
#' @return a list (pvalue=, observed_score=, nb_permutations=)
#' @export
simpleGenotypicSumScoresPvalue <- function(nb_permutations,genotypes_list,  phenotypes,
		covariables=data.frame(), regCoeff=0) {

	n <- length(phenotypes)

	observed <- simpleGenotypicSumScores(genotypes_list, phenotypes, covariables, regCoeff)

	# permute phenotypes and compute corresponding score
	scores <- sapply(1:nb_permutations, function(x) {
				# permute the phenotypes (as a side-effect) by covariables
				by(1:n, covariables, function(x) phenotypes[x] <<- sample(phenotypes[x]) )

				simpleGenotypicSumScores(genotypes_list, phenotypes, covariables, regCoeff)
			})

	# compute pvalue
	nb_better <- sum(scores >= observed)
	pvalue <- (nb_better+1)/(nb_permutations+1)
	return( list(pvalue=pvalue, observed=observed, nb_better=nb_better, nb_permutations=nb_permutations) )

}