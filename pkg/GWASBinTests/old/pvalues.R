#' compute a pvalue using a score function taking a data frame of variables
#'
#' This will repeatedly permute the variables using \code{\link{permuteVariables}}
#' and compute the corresponding score.
#'

#' @param a data frame with at least one variable
#' @param nb_permutation the number of permutations to use
#' @param score the score function
#'
#' @return a list (pvalue=, observed_score=, nb_better=, nb_permutations=)
#'
#' @export
estimatePvalue <- function(variables, nb_permutation, score) {
	observed_score <- score(variables)
	scores <- sapply(1:nb_permutation, function(x) score( permuteVariables(variables) ) )

	nb_better <- sum(scores >= observed_score)
	pvalue <-  	(nb_better+1)/(nb_permutation+1)
	return( list(pvalue=pvalue, observed_score=observed_score, nb_better=nb_better, nb_permutations=nb_permutation) )
}


#' permute the first variable of a data frame keeping the other variables constant
#'
#' Given a data frame with one main variable (the first one)
#' and some co-variables, this function will permute the data frame by preserving
#' the co-variables structure
#
#' @param a data frame with at least one variable
#'
#' @export
permuteVariables <- function(variables) {
	n <- length(variables[,1])
	il <- by(1:n, variables[-1], function(i) i)
	df <- permuteByGroup(variables, il)
	return(df)
}

#' permute lists of subsets of rows
#'
#' Given a list of vector of indices, will permute
#' the rows of x
#'
#' @param x a data structure with rows
#' @param indices_list a list of vector of indices
#'
#' @return the permuted x
#'
#' @export
permuteByGroup <- function(x, indices_list) {
	for (indices in indices_list) {
		if ( ! is.null(indices) )  # indices can be null: correspond to combination of co-variables values not observed
			x[sample(indices),] <- x[indices,]
	}
	return(x)
}

