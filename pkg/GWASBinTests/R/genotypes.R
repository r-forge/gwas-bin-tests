# functions related to genotypes

#' fetch the genotypes from a \strong{GenABEL} dataset in our format
#'
#' @param gws the \strong{GenABEL} dataset
#' @param index the index of the SNP. If NULL (the default), will fetch genotypes for all the SNPs
#'
#' @return the genotypes, an integer matrix with values meaning: \describe{
#' 		\item{homozygous1 (aa) }{0}
#' 		\item{heterozygous (Aa) }{1}
#' 		\item{homozygous2 (AA) }{2}
#' 		\item{missing (\eqn{\emptyset}) }{-1}
#' }
#'
#'
#' @examples
#' data(srdta)
#' # one index
#' g <- fetchGenotypes(srdta,1)
#'
#' # index range
#' gs <- fetchGenotypes(srdta, 1:10)
#'
#' # names
#' gs <- fetchGenotypes(srdta, c("rs150", "rs179"))
#'
#' # all
#' gs <- fetchGenotypes(srdta)
#'
#' @export
fetchGenotypes <- function(gws, index=NULL) {
	if ( is.null(index) )
		g <- as.double(gws[,])
	else
		g <- as.double(gws[,index])
	gi <- as.integer(g)
	dim(gi) <- dim(g)
	return( convertGenotypes(gi) )
}

#' fetch the genotypes from a \strong{GenABEL} dataset in our format as list of vectors
#'
#' See \code{\link{fetchGenotypes}}
#'
#' @param gws the \strong{GenABEL} dataset
#' @param index the index of the SNP(s). If NULL (the default), will fetch genotypes for all the SNPs
#'
#' @return the genotypes, a list of genotypes integer vectors
#' @export
fetchGenotypesAsList <- function(gws, index=NULL) {
	m <- fetchGenotypes(gws, index)
	l <- as.list( as.data.frame(m))
	names(l) <- snpnames(gws)[index]
	return(l)
}

#' convert and check the genotypes to our simple pure-R implementation conventions
#'
#' The genotypes must be an integer vector.
#' The code for missing data is -1, so all NAs will be converted into -1
#'
#' All other values will raise an error.
#'
#' This format of genotypes is expected by all the simple_* functions.
#' To extract genotypes from a Gws, rather use the \code{\link{fetchGenotypes}} and \code{\link{fetchGenotypesAsList}}
#' functions that will call this function anyway.
#'
#' @param genotypes a vector (or matrix) of integers
#' 	with values meaning: \describe{
#' 		\item{homozygous1 (AA) }{0}
#' 		\item{heterozygous (AB) }{1}
#' 		\item{homozygous2 (BB) }{2}
#' 		\item{missing (\eqn{\emptyset}) }{\code{NA} or -1}
#' }. There is no assumption about which allele is minor.
#'
#' @return the cleaned and checked integer vector of genotypes, with values in -1:2
#'
#' @examples
#' data(srdta)
#' 	gs <- convertGenotypes( as.integer( as.double(srdta[,1:20]) ) )
#'
#' @export
convertGenotypes <- function(genotypes) {

	if ( length(genotypes) == 0 )
		stop("got no genotypes")

	if ( ! is.integer(genotypes) )
		stop("genotypes must be an integer vector")

	# N.B: -1L instead of -1 is very important, otherwise it turns the vector to double
	genotypes[is.na(genotypes)] <- -1L

	# check range
	if( any( genotypes < -1 | genotypes > 2 ) )
		stop("Bad genotype values, valid range is -1:2")

	return(genotypes)
}

#' compute the missing rate from the genotypes of a marker
#'
#' The missing rate (\eqn{\beta}) is the proportion of missing data (i.e. NA)
#' among the genotypes
#'
#' @param genotypes a vector of integers: see \code{\link{convertGenotypes}}
#'
#' @return the missing rate
#'
#' @examples
#' data(srdta)
#' beta <- missingRateFromGenotypes( fetchGenotypes(srdta, 1) )
#'
#' @export
missingRateFromGenotypes <- function(genotypes) {
	genotypes <- convertGenotypes(genotypes)

	beta <- sum( genotypes < 0 ) / length(genotypes)
	return(beta)
}
