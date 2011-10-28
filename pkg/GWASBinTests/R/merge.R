#' Merge several genome wide association studies into a single dataset
#'
#' Merge two association studies datasets as \code{\link[GenABEL]{gwaa.data-class}} objects into a single one.
#' The 'gws' column in phenodata keeps track of the sample study of origin.
#'
#' \strong{N.B} the sample ids have to be distinct !!!
#'
#' \strong{N.B} When using intersected_snps_only=FALSE, genotypes for samples for SNPs not in dataset
#' will be set to NA.
#'
#'
#'
#' When doing a GWASBinTests analysis, the method takes into account the study of origin via the \strong{gws}
#' variable by permuting the phenotypes of the labels intra-study.
#'
#' \strong{N.B}: the function is named \code{mergeGws} and not \code{merge.gws} to avoid \code{R CMD check}
#' warnings about \emph{S3 generic/method consistency} and the generic \code{\link{merge}} S3 function.
#'
#' @param gwsA First genome wide association study as a \code{\link[GenABEL]{gwaa.data-class}} object.
#' @param gwsB Second genome wide association study as a \code{\link[GenABEL]{gwaa.data-class}} object.
#' @param intersected_snps_only merge only SNPs shared by gwsA and gwsB
#'
#' @return 	A genome wide association dataset as a \code{\link[GenABEL]{gwaa.data-class}} object.
#' 	Only the following fields are merged:
#'
#' 			\item{pop}{ Case/Control status;}
#' 			\item{sex}{ Sex;}
#' 			\item{gws}{ Genome wide study identifier.}
#'
#' If the \code{gws} field is absent in one of the studies,
#' it is added with a value different from the ones found in the other study.
#'
#' @export
mergeGws <- function( gwsA, gwsB, intersected_snps_only=TRUE)
{
	# check the sample ids are distinct
	if ( any( gwsA@phdata$id %in%  gwsB@phdata$id) )
		stop("ERROR: There are some common sample ids in the 2 gws: ",
				paste(gwsA@phdata$id[gwsA@phdata$id %in%  gwsB@phdata$id], collapse = ", "))

	#It is not a problem to change data as all is done by copy...
	if (is.null(gwsA@phdata$gws))
		gwsA@phdata$gws <- 1L

	maxA <- as.integer(max(gwsA@phdata$gws)) + 1L

	if ( is.null(gwsB@phdata$gws) ) {
		gwsB@phdata$gws <- maxA
	} else {
		gwsB@phdata$gws <- maxA + gwsB@phdata$gws
	}

	gwsAB <- merge.gwaa.data(gwsA,gwsB,intersected_snps_only=intersected_snps_only);
	# the merge should create var.x and var.y for each variable var

	for ( var in c("sex", "pop", "gws") ) {
		varx <- paste(var, ".x", sep="")
		vary <- paste(var, ".y", sep="")

		if ( varx %in% names(gwsAB@phdata) ) {

			varx_na <- is.na(gwsAB@phdata[varx])

			gwsAB@phdata[var] <- gwsAB@phdata[varx]
			gwsAB@phdata[var][varx_na] <- gwsAB@phdata[vary][varx_na]


			# the following completely removes the variables
			gwsAB@phdata[varx] <- NULL
			gwsAB@phdata[vary] <- NULL
		}

	}

	return(gwsAB );
}


