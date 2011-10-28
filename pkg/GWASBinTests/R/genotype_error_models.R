

#' compute the \emph{simple} genotype error model for a marker (SNP)
#'
#' See \link{concept#genotype.error.model}.
#'
#' This model is much simpler than the affymetrix one, in that
#' in does not assume that there is a bias on the heterozygous genotypes.
#'
#' This will generate the following genotype error model:
#'
#' \deqn{
#'  P( \mathcal{O}^j \backslash  \mathcal{G}^j ) =
#'  \left(
#'  \begin{array}{r|r|r|r}
#'  \mathcal{O}^j \backslash  \mathcal{G}^j  & aa & Aa & AA\\
#'  \hline
#'
#'  aa & (1-\beta)*(1-\alpha) & (1-\beta)*\alpha/2 &  (1-\beta)*\alpha/2 \\
#'  \hline
#'  Aa & (1-\beta)*\alpha/2 & (1-\beta)*(1-\alpha) & (1-\beta)*\alpha/2 \\
#'  \hline
#'  AA & 1-\beta)*\alpha/2 & 1-\beta)*\alpha/2 & (1-\beta)*(1-\alpha) \\
#'  \hline
#'  \emptyset &  \beta & \beta & \beta \\
#'  \end{array}
#'  \right)
#' }
#'
#' Some comments: \enumerate{
#' \item \eqn{ (1-\beta) } is the probability for an observed genotype not to be missing.
#' That is why this term is included in all non-missing observed genotype rows.
#' \item \eqn{ (1-\alpha) } is the probability of not doing an error
#' \item the probability of not doing any mistake (the diagonal) for a given genotype is \eqn{(1-\beta)*(1-\alpha)},
#'  which simply means that this is not missing AND there is no error.
#' \item the error is equally split among the two other genotypes
#'
#' }
#'
#' @param alpha the error_rate \eqn{ \alpha }. See \strong{\link{parameters}}
#' @param beta the missing value rate \eqn{ \beta }
#'
#' @return return error model a double
#'
#' @examples
#' m <- simpleGenotypeErrorModel(0.05, 0.01)
#'
#' # error model for first SNP of srdta
#' # by computing missing rate on genotypes
#' data(srdta)
#' beta <- missingRateFromGenotypes( fetchGenotypes(srdta,1) )
#' m <- simpleGenotypeErrorModel(0.05, 0.01)
#'
#' prob <- m["aa","Aa"]
#'
#' # error model for no errors !
#' no_errors <- simpleGenotypeErrorModel(0,0)
#'
#'
#' @export
simpleGenotypeErrorModel <- function(alpha, beta) {
	if ( alpha < 0 || alpha > 1)
		stop("Bad alpha_rate value")
	if ( beta < 0 || beta > 1)
		stop("Bad beta_rate value")

	# prob for having different genotypes
	p_other <- as.double( (1-beta)*alpha/2 )

	m <- matrix(p_other,nrow = 4, ncol = 3)
	rownames(m) <- c("aa", "Aa", "AA", "ZZ")
	colnames(m) <- c("aa", "Aa", "AA")

	# prob for no error
	p_same <- (1-beta)*(1-alpha)
	diag(m) <- p_same

	m["ZZ",] <- beta

	return(m)
}

#