#' Save a genome wide association study dataset (genotypes + phenotypes)
#'
#' Save a genome wide association study dataset in the internal GenABEL format.
#' It is just a proxy for the \code{\link[GenABEL]{save.gwaa.data}} function of \code{\link[GenABEL]{GenABEL}}.
#'
#' The dataset will be stored in two files named \code{basename}\strong{.gag}
#'  and \code{basename}\strong{.gap}
#' for the genotypes and the phenotypes.
#'
#'
#' @param gws A genome wide association dataset as a \code{\link[GenABEL]{gwaa.data-class}} object.
#' @param basename Full path to files without extensions.
#'   The function adds successively \code{.gag} and \code{.gap} to the path in order to have the full paths to files.
#' @param verbose if TRUE, suppress the output of the GenABEL \code{\link[GenABEL]{save.gwaa.data}} function.
#'
#' @return No value returned
#'
#' @export
#'

saveGws <- function(gws, basename, verbose = FALSE) {
	save <- function() save.gwaa.data(gws, phenofile = paste(basename,".gap", sep=""),
				genofile = paste(basename,".gag", sep=""), human = FALSE)

	if ( verbose )
		gws <- save()
	else
		output <- capture.output( save() )

	return(NULL)
}