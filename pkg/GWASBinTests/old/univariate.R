

#' compute the contingency table for one SNP and some variables
#'
#' The NAs in the variables are not taken into account, i.e. the corresponding
#' samples are nout counted in the table
#'
#' @param genotypes a vector of integers: see \code{\link{convertGenotypes}}
#' @param vars a data frame of variables.s
#'
#' @return the table
#'
#'
#' @export
univariateTable <- function(genotypes, vars=NULL) {
	if ( ! is.null(vars) && length(genotypes) != length(vars[,1]))
		stop("Different length between genotypes and phenotypic data")

	f <- factor(genotypes, levels=c(0:2, -1), labels=c("aa", "Aa", "AA", "ZZ"))

	if ( ! is.null(vars) )
		df <- data.frame(f, vars, check.rows = TRUE, check.names = TRUE)
	else
		df <- data.frame(f, check.rows = TRUE, check.names = TRUE)
	t <- table(df)

	return(t)
}

#' compute the genotype probabilities
#'
#' This function will compute the estimate \eqn{ \widehat{P}(\mathcal{S}, (\mathcal{V})_m, \mathcal{G}) }
#' of the distribution of \eqn{ P(\mathcal{S}, (\mathcal{V})_m, \mathcal{G}) } for a given marker.
#'
#' Each genotype is a realization of this distribution, so that we can compute the estimation from the contingency tables:
#' \deqn{
#'  	\widehat{P}(\mathcal{S}, (\mathcal{V})_m, \mathcal{G}) =
#'  \frac{ n(\mathcal{S},(\mathcal{V})_m,\mathcal{G}) }{ n(\mathcal{S},(\mathcal{V})_m,\oplus) }
#' }
#'
#' But the \eqn{\mathcal{G}} are hidden variables, so the probabilities will be estimated
#' from the observed genotype \eqn{\mathcal{O}} frequencies through the removal of patients
#'  with missing genotypes:

#'
#' \deqn{
#'  	\widehat{P}(\mathcal{S}, (\mathcal{V})_m, \mathcal{G}) = \widehat{P}(\mathcal{S}, (\mathcal{V})_m, \mathcal{O})
#'  = \frac{ n(\mathcal{S},(\mathcal{V})_m,\mathcal{O}) + C }{ n(\mathcal{S},(\mathcal{V})_m, \mathcal{O} \ne \emptyset) + 3 \times C}
#' }
#' where \eqn{C=r . \bar{n}}  is the regularization constant
#' (see \strong{\link{concept#regularization.of.table.estimates}} )
#'
#' In other words, the genotype frequencies are estimated from the contingency table of observed genotypes
#' without the cells corresponding to missing data \eqn{\emptyset}, and with a constant added to each cell if the regularization
#' coefficient \strong{r} is not null.
#'
#'
#' For example, in the case of a binary phenotype (S=0 or 1) and no co-variable, with the following notations:
#'
#' \eqn{\begin{tabular}{|c||c|c|}
#' \hline
#' O & S=0 & S=1\tabularnewline
#' \hline
#' \hline
#' aa & a & e\tabularnewline
#' \hline
#' Aa & b & e\tabularnewline
#' \hline
#' AA & c & f\tabularnewline
#' \hline
#' $\emptyset $ & x & y\tabularnewline
#' \hline
#' \end{tabular}
#' }
#'
#' we have:
#' \deqn{
#'  	\widehat{P}(\mathcal{S}=0, \mathcal{G}=aa)
#'  = \frac{ a + C }{ a + b + c + 3 \times C }
#' }
#' with \eqn{ C=r . \frac{a+b+c}{3}}
#'
#' @param table an object of class table, with at least 2 dimensions (genotypes and phenotype \eqn{ \mathcal{S} }).
#' See \code{\link{univariateTable}}
#' @param   regCoeff 	regularization rate: will be used to compute the constant to add
#' 	in each cell of the contingency table to regularize estimators. This constant will
#' 	be = rate*total/nb_of_cells. See \strong{\link{concept#regularization.of.table.estimates}}.
#'
#' @return an array of doubles of same dimensions than the table, i.e \eqn{\widehat{P}(\mathcal{S}, (\mathcal{V})_m, \mathcal{G})}
#' N.B for convenience, the array will have a first dimension of length 4, even if there are only 3 genotypes. The 4th values
#' corresponding to the O=ZZ will be filled with NA.
#' If there is no genotype for a category then the probabilty can not be calculated and so the probability is NaN
#'
#' @export
univariateGenotypeProbabilities<- function(table, regCoeff=0) {
	d <- length(dim(table))


	# instead of removing the ZZ row table["ZZ",..] I will put NA
	table_nozz <- table
	zz_indices <- slice.index(table_nozz,1) == 4
	table_nozz[ zz_indices ] <- NA

#	# remove the part corresponding to the missing data phenotype
#	table_nozz <- table[slice.index(table,1) < 4] # 4 is ZZ
#	# fix dimensions and dimension names
#	dims <- dim(table)
#	dims[1] <- 3
#	dnn <- dimnames(table)
#	dnn[[1]] <- dnn[[1]][-4]
#	dim(table_nozz) <- dims
#	dimnames(table_nozz) <- dnn

	if ( regCoeff != 0) {
		sums <- apply(table_nozz, 2:d, sum, na.rm = TRUE)
		reg_constants <- sums*regCoeff/3
		table_nozz <- sweep(table_nozz, 2:d, reg_constants, FUN="+")
	}

	sums <- apply(table_nozz, 2:d, sum, na.rm = TRUE)
	res <- sweep(table_nozz, 2:d, sums, FUN = "/")

	# if we have some zeroes in sums, we have introduced some NaN
	# but they will not be used

#	# fill the row of ZZ with NAs because theses values should not be used
#	res[ zz_indices ] <- NA

	return(unclass(res))
}

#' compute the univariate log likelihood for a given marker
#'
#'
#' See \strong{\link{concept#naive.likelihood}}
#'
#' This function computes the log of the naive L3 likelihood of a marker j.
#' \deqn{
#' 	L3(Marker_j) = \prod_{all \ patients \ i} \ \sum_{g \in \{aa, Aa, AA\} } p( o^j(i) | g) . p(g | (v(i))_m
#' }
#'
#' we can see that it only depends on \eqn{ o^j(i) } and \eqn{ (v(i))_m }
#' so that we can write:
#'
#' \deqn{
#' 	L3(Marker_j) = \prod_{ (o, (v)_m) \in \{aa, Aa, AA,\emptyset\} \times \{(V)_m\} }
#' 							[ \sum_{g \in {aa, Aa, AA}} p( o| g) . p(g | (v)_m) ) \; ] ^{ n( (v)_m, o) }
#' }
#'
#' Hence:
#'
#' \deqn{
#' 	log(L3(Marker_j)) = \sum_{ (o, (v)_m) \in \{aa, Aa, AA,\emptyset\} \times \{(V)_m\} \backslash n( (v)_m, o) \ne 0 }
#' 							n( (v)_m, o) . log( \sum_{g \in {aa, Aa, AA}} p( o| g) . p(g | (v)_m) ) )
#' }
#'
#' In order to compute this log-likelihood, we need: \enumerate{
#' \item the contingency table of observed genotypes and a set of patient variables (V)m that gives us the \eqn{ n( (v)_m, o) }
#' \item  an array of genotype probabilities on the same variables: \eqn{  p(g | (v)_m) ) }
#' \item a genotyping error model: \eqn{P( O | G)}
#'
#' }
#'
#' @param gct the genotype contingency table, see \code{\link{univariateTable}}
#' @param probs an array of genotype probabilities, see \code{\link{univariateGenotypeProbabilities}}
#' @param error_model a genotyping error model, see \link{concept#genotype.error.model}, \code{\link{simpleGenotypeErrorModel}}
#'
#' @return the log likelihood numeric value. If there is no genotype for a category, the likelihood will not be defined
#' 	and so it will return NaN
#' @export
univariateLogLikelihoodForMarker <- function(gct, probs, error_model) {
	d <- length(dim(probs))

	# compute the probabilities of each observed genotype
	observed_genotype_probs <- apply(probs, 2:d, function(x) error_model %*% x[-4] )

	# switch to log
	log_obs_genotype_probs <- log(observed_genotype_probs)

	# make the term-by-term product, then sum over the non-empty cells of the table
	non_zero <- which(gct > 0) # non-empty cells
	llh <- sum( (gct * log_obs_genotype_probs)[non_zero] )
	return(llh)
}



#' compute the univariate likelihood Ratio for a dataset/bin
#'
#' See \strong{\link{concept#likelihood.ratio}} and \strong{\link{concept#naive.likelihood}}.
#' For hypotheses H1 and H0, it computes the sum of log-likelihoods for all markers
#' (using \code{\link{univariateLogLikelihoodForMarker}})
#' then substract the sums.
#'
#' More precisely:
#' \deqn{
#' 	logL3Ratio( Bin_b ) =  \sum_{all \ markers \ j} \  log(L3_{H_1}( Marker_j))
#' 							- \sum_{all \ markers \ j} \  log(L3_{H_0}( Marker_j))
#' }
#'
#' @param genotypes a matrix of integers: see \code{\link{fetchGenotypes}}
#' @param variables a data frame of phenotype and optional co-variables: \eqn{ (\mathcal{S},(\mathcal{V})_m) }
#' 	the phenotype must be the first variable of the data frame
#' @param error_model a genotyping error model, see \link{concept#genotype error model}, \code{\link{simpleGenotypeErrorModel}}

#' @param  ... will be sent to \code{\link{univariateGenotypeProbabilities}} optional arguments
#'
#' @return the log-likelihood ratio
#' @export
univariateLikelihoodRatio <- function(genotypes, variables, error_model, ...) {
	if ( ! is.data.frame(variables) )
		stop("variables must be a data frame")
	if ( length( variables ) < 1 )
		stop("variables data frame must contain at least oen variable: the phenotype")

	h1 <- variables
	h0 <- h1
	h0[1] <- 0 # nulls out the phenotypes

	hypothesis <- list(h1=h1,h0=h0)

	all_llr <- apply(genotypes, 2, function(genos) {
		llhs <- lapply(hypothesis, function(h) {
					ct <- univariateTable(genos, h)
					probs <- univariateGenotypeProbabilities(ct, ...)
					llh <- univariateLogLikelihoodForMarker(ct, probs, error_model)
				})

		return(llhs[[1]] - llhs[[2]])
	})

	llr <- sum(all_llr)

	return(llr)
}

#' compute the univariate likelihood Ratio pvalue for a dataset/bin using permutations
#'
#' See \code{\link{univariateLikelihoodRatio}}
#'
#' @param nb_permutations the number of permutations to do
#' @param genotypes a matrix of integers: see \code{\link{fetchGenotypes}}
#' @param variables a data frame of phenotype and optional co-variables: \eqn{ (\mathcal{S},(\mathcal{V})_m) }
#' 	the phenotype must be the first variable of the data frame
#' @param error_model a genotyping error model, see \link{concept#genotype error model}, \code{\link{simpleGenotypeErrorModel}}

#' @param  ... will be sent to \code{\link{univariateGenotypeProbabilities}} optional arguments
#'
#' @return a list (pvalue=, observed_score=, nb_permutations=)
#' @export
univariatePvalue <- function(nb_permutations, genotypes, variables, error_model, ...) {
	score <- function(vars) {
		univariateLikelihoodRatio(genotypes, vars, error_model, ...)
	}
	res <- estimatePvalue(variables, nb_permutations, score)
	return(res)
}
