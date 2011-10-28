# utilities to test/compare the different implementations

.generic_compare_vs_nnbc_pvalues <- function(gws, nnbc_params, phenotypes, covariables, bins,
		nnbc_type, score_function, simple_permutations, pvalue_function) {


	df <- .generic_compare_vs_nnbc(gws, nnbc_params, phenotypes, covariables, bins,
			nnbc_type, score_function)

	# compute pvalues using the bins from nnbc results
	if ( nrow(df) == 0)
		return(data.frame())

	if ( simple_permutations > 0 ) {

		pvs_simple <- apply(df, 1, function(b) {
					gws_bin <- gwsForBin(gws, b['chr'], b['start'], b['end'])
					genotypes <- fetchGenotypesAsList(gws_bin)
					pv <- pvalue_function(simple_permutations, genotypes, phenotypes, covariables)
					return(pv$pvalue)
				})

		df[paste("simple_", nnbc_type, "_pv", sep="")] <- pvs_simple

		# compute a confidence interval

		df$in_conf_int <- mapply(function(pv, pve) {
					ci <- binconf(pv*simple_permutations, simple_permutations, 0.001)

					return(pve >= ci[2] && pve <= ci[3])
				},
				df[,paste("pv_", nnbc_type, sep="")], pvs_simple)

	}

	return(df)
}

.generic_compare_vs_nnbc <- function(gws, nnbc_params, phenotypes, covariables, bins,
		nnbc_type, score_function) {
	if (nnbc_params$use_affymetrix_model)
		stop("use_affymetrix_model not supported in compareUnivariateLikelihood")

	nnbc_params$types <-  nnbc_type

	fake_bin <- FALSE
	if ( is.null(bins)) {
		fake_bin <- TRUE
		# make a bin fitting exactly gws
		chrs <- unique(chromosome(gws))
		if ( length(chrs) > 1 )
			stop("all snps must be on the same chromosome")
		min_max <- range(map(gws))
		# N.B: the order is paramount !!!!!!!
		bins <- data.frame(chr_name=chrs, start=min_max[1], end=min_max[2],bin_index=c(1L) )
	}

	# run nnbc (c++)

	# we want to force the nnbc c++ to use the exact same covariables
	# so we null out the gws column if any
	gws@phdata$pop <- phenotypes
	if ( "gws" %in% names(gws@phdata) )
		gws@phdata$gws <- rep(0, nrow(gws@phdata))
#	nnbc_params$nb_permutations <- 0

	res_nnbc <- processGws(gws, bins, nnbc_params, covariables=covariables)

	if (fake_bin) {
		genotypes <- fetchGenotypesAsList(gws)
		res_simple <- score_function(genotypes, phenotypes, covariables)
	} else {
		# use the bins from nnbc results
		res_nnbc <- res_nnbc[res_nnbc$nb_snps > 1,]
		if ( nrow(res_nnbc) == 0)
			return(data.frame())

		res_simple <- apply(res_nnbc, 1, function(b) {
					gws_bin <- gwsForBin(gws, b['chr'], b['start'], b['end'])
					genotypes <- fetchGenotypesAsList(gws_bin)
					return(score_function(genotypes, phenotypes, covariables))
				})
	}

	res_nnbc[paste("simple_", nnbc_type, "_score", sep="")] <- res_simple
	return(res_nnbc)
}

#' compute the univariateLikelihood using both implementations
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param nnbc_params the NNBC parameters, see \code{\link{parameters}}
#' @param phenotypes an optional vector of phenotypes
#' @param bins an optional dataframe of bins. if not given it will use a single bin
#'        containing all the dataset
#' @return a data frame with columns: bin chr start end nb_snps nb_permutations pv_divariate score_divariate fdr_divariate simple_divariate_llh
.compareUnivariateLikelihood <- function(gws, nnbc_params, phenotypes=gws@phdata$pop, covariables=data.frame(), bins=NULL) {
	.compareUnivariateLikelihoodPvalues(gws, nnbc_params, 0, phenotypes, covariables, bins)
}


#' compute the H1 univariateLikelihood pvalues using both implementations
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param nnbc_params the NNBC parameters, see \code{\link{parameters}}
#' @param simple_permutations the number of permutations for the R simple implementation
#' @param phenotypes an optional vector of phenotypes
#' @param bins an optional dataframe of bins. if not given it will use a single bin
#'        containing all the dataset
#' @return a list of results
.compareUnivariateLikelihoodPvalues <- function(gws, nnbc_params, simple_permutations,
		phenotypes=gws@phdata$pop, covariables=data.frame(), bins=NULL ) {

	univariate_score <- function(genotypes, phenotypes, covariables) {
		if ( length(covariables) )
			variables <- data.frame(phenotypes, covariables)
		else
			variables <- data.frame(phenotypes)
		error_models <- lapply(genotypes, function(g) simpleGenotypeErrorModel(nnbc_params$max_error_rate, missingRateFromGenotypes(g)))
		simpleUnivariateLogLikelihood(genotypes, variables, error_models)
	}

	univariate_pvalue <- function(simple_permutations, genotypes, phenotypes, covariables) {
		error_models <- lapply(genotypes, function(g) simpleGenotypeErrorModel(nnbc_params$max_error_rate, missingRateFromGenotypes(g)))
		simpleUnivariatePvalue(simple_permutations, genotypes, phenotypes, error_models, covariables)
	}

	.generic_compare_vs_nnbc_pvalues(gws, nnbc_params, phenotypes, covariables, bins, "univariate", univariate_score,
			simple_permutations, univariate_pvalue)

}


#' compute the H1 DivariateLikelihood pvalues using both implementations
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param nnbc_params the NNBC parameters, see \code{\link{parameters}}
#' @param simple_permutations the number of permutations for the R simple implementation
#' @param phenotypes an optional vector of phenotypes
#' @param bins an optional dataframe of bins. if not given it will use a single bin
#'        containing all the dataset
#' @return a list of results

.compareDivariateLikelihoodPvalues <- function(gws, nnbc_params, simple_permutations,
		phenotypes=gws@phdata$pop, covariables=data.frame(), bins=NULL) {


	divariate_score <- function(genotypes, phenotypes, covariables) {
		if ( length(covariables) )
			variables <- data.frame(phenotypes, covariables)
		else
			variables <- data.frame(phenotypes)
		error_models <- lapply(genotypes, function(g) simpleGenotypeErrorModel(nnbc_params$max_error_rate, missingRateFromGenotypes(g)))
		simpleDivariateLogLikelihood(genotypes, variables, error_models)
	}

	divariate_pvalue <- function(simple_permutations, genotypes, phenotypes, covariables) {
		error_models <- lapply(genotypes, function(g) simpleGenotypeErrorModel(nnbc_params$max_error_rate, missingRateFromGenotypes(g)))
		simpleDivariatePvalue(simple_permutations, genotypes, phenotypes, error_models, covariables)
	}

	.generic_compare_vs_nnbc_pvalues(gws, nnbc_params, phenotypes, covariables, bins, "divariate", divariate_score,
			simple_permutations, divariate_pvalue)

}



#' compute the H1 DivariateLikelihood  using both implementations
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param nnbc_params the NNBC parameters, see \code{\link{parameters}}
#' @param phenotypes an optional vector of phenotypes
#' @param bins an optional dataframe of bins. if not given it will use a single bin
#'        containing all the dataset
#' @return a data frame with columns: bin chr start end nb_snps nb_permutations pv_divariate score_divariate fdr_divariate simple_divariate_llh

.compareDivariateLikelihood <- function(gws, nnbc_params,
		phenotypes=gws@phdata$pop, covariables=data.frame(), bins=NULL ) {

	.compareDivariateLikelihoodPvalues(gws, nnbc_params, 0, phenotypes, covariables, bins)

}

#' compute the allelic sum score using both implementations
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param nnbc_params the NNBC parameters, see \code{\link{parameters}}
#' @param phenotypes an optional vector of phenotypes
#' @param bins an optional dataframe of bins. if not given it will use a single bin
#'        containing all the dataset
#' @return a data frame with columnd: bin chr start end nb_snps nb_permutations pv_divariate score_divariate fdr_divariate simple_divariate_llh
.compareAllelicSumScores <- function(gws, nnbc_params, phenotypes=gws@phdata$pop, covariables=data.frame(), bins=NULL) {
		.compareAllelicSumScoresPvalues(gws, nnbc_params, 0, phenotypes, covariables, bins)
}

#' compute the allelic sum score using both implementations
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param nnbc_params the NNBC parameters, see \code{\link{parameters}}
#' @param phenotypes an optional vector of phenotypes
#' @param bins an optional dataframe of bins. if not given it will use a single bin
#'        containing all the dataset
#' @return a data frame with columnd: bin chr start end nb_snps nb_permutations pv_divariate score_divariate fdr_divariate simple_divariate_llh

.compareAllelicSumScoresPvalues <- function(gws, nnbc_params, simple_permutations,
		phenotypes=gws@phdata$pop, covariables=data.frame(), bins=NULL) {

	.generic_compare_vs_nnbc_pvalues(gws, nnbc_params, phenotypes, covariables, bins, "allelic",
			simpleAllelicSumScores,
			simple_permutations, simpleAllelicSumScoresPvalue)
}

#' compute the genotypic sum score using both implementations
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param nnbc_params the NNBC parameters, see \code{\link{parameters}}
#' @param phenotypes an optional vector of phenotypes
#' @param bins an optional dataframe of bins. if not given it will use a single bin
#'        containing all the dataset
#' @return a data frame with columnd: bin chr start end nb_snps nb_permutations pv_divariate score_divariate fdr_divariate simple_divariate_llh
.compareGenotypicSumScores <- function(gws, nnbc_params,
		phenotypes=gws@phdata$pop, covariables=data.frame(), bins=NULL ) {

	.compareGenotypicSumScoresPvalues(gws, nnbc_params, 0, phenotypes, covariables, bins)
}

#' compute the genotypic sum score using both implementations
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param nnbc_params the NNBC parameters, see \code{\link{parameters}}
#' @param phenotypes an optional vector of phenotypes
#' @param bins an optional dataframe of bins. if not given it will use a single bin
#'        containing all the dataset
#' @return a data frame with columnd: bin chr start end nb_snps nb_permutations pv_divariate score_divariate fdr_divariate simple_divariate_llh

.compareGenotypicSumScoresPvalues <- function(gws, nnbc_params, simple_permutations,
		phenotypes=gws@phdata$pop, covariables=data.frame(), bins=NULL) {

	.generic_compare_vs_nnbc_pvalues(gws, nnbc_params, phenotypes, covariables, bins, "genotypic",
			simpleGenotypicSumScores,
			simple_permutations, simpleGenotypicSumScoresPvalue)
}
