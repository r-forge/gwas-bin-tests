# functions related to bins

#' create a subset of a gws for a given bin
#'
#' @param gws the dataset as a \code{\link[GenABEL]{gwaa.data-class}} object
#' @param chr the chromosome of the bin
#' @param start the start of the bin
#' @param end the end of the bin
#'
#' @return a (sub) dataset as a \code{\link[GenABEL]{gwaa.data-class}} object, or NULL
#' @export
gwsForBin <- function(gws, chr, start, end) {

	gt <- gtdata(gws)
	in_bin <- gt@chromosome == chr & gt@map >= start & gt@map <= end
	#cat(chr,start, end, sum(in_bin), "\n")
	if ( any(in_bin))
		return(gws[,in_bin])
	else
		return(NULL)
}

#' make a bin for every SNP of the dataset
#'
#' @param gws the dataset as a \code{\link[GenABEL]{gwaa.data-class}} object
#'
#' @return a sorted dataframe of bins (chr_name, start, end, bin_index)
#' @export
makeBinsForSnps <- function(gws) {

	df <- data.frame(chr_name=chromosome(gws), start=map(gws), end=map(gws))


	chrs <- convertChromosomeNames(df$chr_name)
	b.o <- order(chrs, df$start)

	df <- df[b.o,]
	df$bin_index <- 1:nrow(df)
	return(df)
}

