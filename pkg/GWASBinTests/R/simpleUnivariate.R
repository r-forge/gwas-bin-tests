

#' compute the univariate likelihood Ratio pvalue for a dataset/bin using permutations
#'
#' See \code{\link{simpleUnivariateLogLikelihoodRatio}}
#' This is a reference implementation, very naive on purpose, optimized for clarity.
#'
#' @param nb_permutations the number of permutations to do
#' @param genotypes_list a list of vector of integers: see \code{\link{fetchGenotypes}}
#' @param phenotypes a vector of phenotypes of same length as the elements of genotypes_list
#' @param error_models list of genotyping error model for each marker, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}
#' if it is not a list, it must be a matrix and this error model will be used for all markers in the dataset
#' @param covariables a data frame of cvariables: \eqn{ (\mathcal{V})_m }
#' @param  regCoeff see \code{\link{simpleUnivariateLogLikelihoodForMarker}}
#'
#' @return a list (pvalue=, observed_score=, nb_permutations=)
#' @export
simpleUnivariatePvalue <- function(nb_permutations, genotypes_list, phenotypes,  error_models, covariables=data.frame(),regCoeff=0) {

	if ( ! is.vector(phenotypes) )
		stop("phenotypes nmust be a vector")
	h1 <- NULL
	if ( length(covariables) )
		h1 <- data.frame(phenotypes=phenotypes, covariables)
	else
		h1 <- data.frame(phenotypes=phenotypes)
	observed <- simpleUnivariateLogLikelihood(genotypes_list, h1,  error_models, regCoeff)


	n <- length(phenotypes)
	scores <- sapply(1:nb_permutations, function(x) {
				# permute the phenotypes (as a side-effect) by covariables
				by(1:n, covariables, function(x) h1$phenotypes[x] <<- sample(h1$phenotypes[x]) )
				simpleUnivariateLogLikelihood(genotypes_list, h1,  error_models, regCoeff)
	})
	nb_better <- sum(scores >= observed)
	pvalue <- (nb_better+1)/(nb_permutations+1)
	return( list(pvalue=pvalue, observed=observed, nb_better=nb_better, nb_permutations=nb_permutations) )
}



#' compute the univariate log likelihood Ratio for a dataset/bin
#'
#' @param genotypes_list a list of vector of integers: see \code{\link{fetchGenotypes}}
#' @param phenotypes a vector of phenotypes of same length as the elements of genotypes_list
#' @param error_models list of genotyping error model for each marker, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}
#' if it is not a list, it must be a matrix and this error model will be used for all markers in the dataset
#' @param covariables a data frame of cvariables: \eqn{ (\mathcal{V})_m }
#' @param  regCoeff see \code{\link{simpleUnivariateLogLikelihoodForMarker}}

#' @return the log-likelihood ratio
#' @export
simpleUnivariateLogLikelihoodRatio <- function(genotypes_list, phenotypes,  error_models, covariables=data.frame(),regCoeff=0) {
	h1 <- NULL
	if ( length(covariables) )
		h1 <- data.frame(phenotypes=phenotypes, covariables)
	else
		h1 <- data.frame(phenotypes=phenotypes)

	llh_h1 <- simpleUnivariateLogLikelihood(genotypes_list, h1, error_models, regCoeff)
	llh_h0 <- simpleUnivariateLogLikelihood(genotypes_list, covariables, error_models, regCoeff)
	return(llh_h1 - llh_h0)
}


#' compute the univariate log likelihood of a dataset
#'
#' This is a reference implementation, very naive on purpose, optimized for clarity.
#'
#' @param genotypes_list a list of vector of integers: see \code{\link{fetchGenotypes}}
#' @param variables a data frame of variables, possibly empty
#' @param error_models list of genotyping error model for each marker, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}
#' if it is not a list, it must be a matrix and this error model will be used for all markers in the dataset
#'
#' @param  regCoeff see \code{\link{simpleUnivariateLogLikelihoodForMarker}}
#'
#' @return the log-likelihood for the set of markers
#' @export
simpleUnivariateLogLikelihood <- function(genotypes_list, variables, error_models, regCoeff=0) {

	nb_markers <- length(genotypes_list)

	if ( ! is.data.frame(variables) )
		stop("variables must be a data frame")

	if ( nb_markers < 1)
		stop("there must be at least genotypes for one marker")

	if ( is.list(error_models)) {
		if ( length(error_models) != nb_markers)
			stop("error_models must have the same length than genotypes_list")
	} else {
		error_models <- rep(list(error_models), nb_markers)
	}
	if ( ! all( sapply(error_models, function(m) is.matrix(m) && dim(m) == c(4,3)) ) )
		stop("bad error models, must be matrices as returned by simpleGenotypeErrorModel")

	sum( mapply( function(g, e) { simpleUnivariateLogLikelihoodForMarker(g, variables, e, regCoeff ) }
					, genotypes_list, error_models) )
}

#' compute the univariate log likelihood for a given marker
#'
#' This is the a very simple implementation, self-contained and easy to check
#' but not performant. It is used as illustration or to check other implementations.
#'
#' @param genotypes a vector of integers: see \code{\link{convertGenotypes}}
#' @param variables a data frame of variables, possibly empty
#' @param error_model a genotyping error model, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}
#' @param  regCoeff 	regularization rate: will be used to compute the constant to add
#' 	in each cell of the contingency table to regularize estimators. This constant will
#' 	be = rate*total/nb_of_cells. See \strong{\link{concept#regularization.of.table.estimates}}.
#'
#' @return the log likelihood numeric value. If there is no genotype for a category, the likelihood will not be defined
#' 	and so it will return NaN
#' @export
simpleUnivariateLogLikelihoodForMarker <- function(genotypes, variables, error_model, regCoeff=0) {

	if ( !is.vector(genotypes) || !is.integer(genotypes))
		stop("genotypes must be a vector of integers !")

	llh3_by_category <- by(genotypes, variables, function(x) {

		###########################################################
		# This is the CORE of the univariate likelihood calculus  #
		###########################################################

		# from here, there is no more variable of co-variable
		# we have only one data: x a vector of observed genotypes for this margin
		# and one error_model

		##  === compute the contingency table  ===
		f <- factor(x, levels=c(0:2, -1), labels=c("aa", "Aa", "AA", "ZZ"))
		ct <- table(f)

		## === compute the estimates of genotype probabilities ===

		ct_no_zz <- ct[-4] # we do not take into account the "ZZ"

		# only missing data ==> we return a llh of zero so that it has no impact on the sum
		if ( sum(ct_no_zz) == 0 )
			return(0)

		# regularization constant if any
		if ( regCoeff != 0) {
			const <- regCoeff * sum(ct_no_zz) / 3 # compute the constant
			ct_no_zz <- ct_no_zz + const # add it to each cell
		}
		probs <- ct_no_zz/ sum(ct_no_zz) # compute the estimates of genotype probabilities

		## == compute the marginal log-likelihood for marker

		# compute the observed genotypes probabilities
		probs_by_observed_genotype <- error_model %*% probs
		dim(probs_by_observed_genotype) <- NULL # trick to get rid of dim(1,n)

		# sum the ni*log(p(oi)) for ni!=0
		non_zero <- which( ct > 0 )

		llh3 <- sum ( (ct * log( probs_by_observed_genotype))[non_zero] )

		# that's all !
		return(llh3)
	})

	return( sum(llh3_by_category) )
}