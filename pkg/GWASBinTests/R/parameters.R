#' make a set of parameters for GWASBinTests
#'
#' It may be used to set some parameter values, and use the defaults values for the others.
#' Moreover some checks are performed and some internal settings computed afterwards.
#' So one should never change a set of parameters after its construction by this function.
#'
#' \subsection{parameters for heuristics}{
#'
#' There are two heuristics implemented, that share the parameter \code{confidence} :
#' \describe{
#' 	\item{the \strong{min_pvalue} heuristic}{
#' 		The goal of this strategy is to avoid speeding permutations on getting a good accuracy on p-values that
#' 		are of no interest. Let's say for instance we are computing p-values for \eqn{10^6} SNPs, we for sure
#' 		are not interested in p-values > 0.1 .
#' 		The estimated p-value follows a binomial distribution, so that we can compute its lower bound for a
#' 		given probability, controlled by the parameter \code{confidence}.
#' 		So if the lower bound of the p-value is greater than \code{min_pvalue} we stop the computation,
#' 		and report the current estimation of the p-value along with the actual number of permutations used to compute it.
#'  }
#'
#' 	\item{the "max_relative_error" heuristic}{
#'	Basically the principle of this heuristic is that the lower the p-value, the higher the number of permutations
#'	it needs to get a good accuracy, and most often we are greatly interested in those small p-values.
#'	What precision do we need ? We could set a maximum value on the radius of the confidence interval
#'	of the estimated p-value, e.g. \eqn{10^6}. But what if the real p-value is 0.01 ?
#'	Then this amount of precision is a waste. And what if the real p-value is \eqn{10^-7} ?
#'	 Then knowing that its confidence interval radius is < \eqn{10^-6} is of not great use.
#'
#'	 So instead the idea is to control the relative error on the approximated p-value
#'  Given the parameters \strong{max_relative_error} and \strong{confidence},
#' 		the computation will stop when the relative error at a thegiven confidence level is lower
#'  than \strong{max_relative_error}
#'
#'   }
#' }
#' }
#'
#'
#' @param types the list of pvalue types to compute. Currently the possible values are
#' 		"univariate", "divariate",  "allelic", "genotypic" (default=all)
#' @param nb_permutations the (maximum) number of permutations to use to compute the pvalues.
#' @param verbosity the amount of verbosity: 0=none, 1=verbose
#' @param seed the seed for the NNBC internal Random Number Generator. Useful to reproduce results.
#' @param max_error_rate The maximum genotyping error rate: \eqn{\alpha}
#' For a given marker, the \strong{error rate} is the probability
#' of having an incorrect (but not missing) observed genotype,
#' i.e. \eqn{ P(\mathcal{O} \ne \mathcal{G} | \mathcal{O} \ne \emptyset) }
#'
#' Hence the maximum error rate \eqn{\alpha} 	is a higher bound of the error rates of all markers,
#' i.e.
#' \deqn{ \forall j, P(\mathcal{O}^j \ne \mathcal{G}^j | \mathcal{O}^j \ne \emptyset) \le \alpha }
#'
#' This \strong{maximum error rate}, sometimes called simply \strong{error rate} is estimated during external
#'  comparison of genotyping technologies. We often use 0.05 as a default value.#'
#'
#' @param use_affymetrix_model Use the \emph{Affymetrix} genotyping error model.
#' @param min_pvalue The minimum "interesting" pvalue. Pvalues above this threshold (with confidence \code{confidence})
#'		may be computed with less permutations.
#' @param confidence value that controls the threshold on probability that a given pvalue computed with a given
#' 		number of permutations is above the \code{min_pvalue} threshold, and so do not need to be computed with
#' 		more permutations
#' @param max_relative_error stop the permutations as soon as the relative error on the palues is below
#'        that threshold at \code{confidence} level. More info in the \emph{Details} section below.
#'
#' @param regularizeEstimators Add a regularization constant to contingency tables.
#' @param excludeX If set, X chromosome is excluded from the study.
#' @param excludeMalesonXChr If set, male patients are excluded for analysis of bins on the X chromosome.
#' @param threads the number of threads to use (if \strong{OPENMP} is supported) to speed up the computation
#'  (defaults to the number of cpus/cores detected on the system).
#'
#' @return a named list with all parameters defined to their default values
#'
#' @examples
#' 	params <- parameters(
#'	  seed = 0,
#'	  verbosity = 0,
#'	  nb_permutations = 100)
#'
#'
#' @export
parameters <- function(
		types = NULL,
		nb_permutations=1000,
		verbosity=0,
		seed=as.integer(Sys.time()),
		max_error_rate=0.05,
		use_affymetrix_model=FALSE,
		min_pvalue = -1,
		confidence = 0.999,
		max_relative_error = 0,
		regularizeEstimators = FALSE,
		excludeX = FALSE,
		excludeMalesonXChr = FALSE,
		threads = 0) {

	params <- list(
			types = types,
			nb_permutations = nb_permutations,
			verbosity = verbosity,
			seed = seed,
			max_error_rate = max_error_rate,
			use_affymetrix_model = use_affymetrix_model,
			min_pvalue = min_pvalue,
			confidence = confidence,
			max_relative_error = max_relative_error,
			regularizeEstimators = regularizeEstimators,
			excludeX = excludeX,
			excludeMalesonXChr = excludeMalesonXChr,
			threads = threads
	)

	return(.checkParameters(params))
}



# ' Check the NNBC parameters set.
# '
# ' Mainly for internal use. May be used to check the integrity of the parameter set
# ' after some changes. In that case one should use the returned value as the
# ' correct parameter set.
# '
# ' @param params the parameters, as returned by \code{\link{parameters}}
# ' @return the possibly corrected set of parameters
# '

.checkParameters <- function(params) {

	TYPES <- c("univariate","divariate", "allelic", "genotypic")

	if(! is.null(params$types) && !all(params$types %in% TYPES) )
		stop("Unknown type(s):", paste(params$types[!params$types %in% TYPES], collapse=", "))

	as_int_params <- c("use_affymetrix_model", "regularizeEstimators", "excludeX", "excludeMalesonXChr",
			 "verbosity", "seed", "threads" )
	params[as_int_params] <- as.integer( params[as_int_params] )

	params$types_as_string <-  paste(params$types, collapse=" ")
	return(params)
}

