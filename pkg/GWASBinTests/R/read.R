

# check the phenodata dataframe of a GenABEL dataset
.checkPhenoData <- function(phdata) {
	# mandatory variables
	names <- names( phdata )
	for ( var in c("id", "sex", "pop") ) {
		if ( ! var %in% names )
			stop("mandatory variable ", var, " not found in phdata !")
	}

	# check integer variables
	for ( var in c("sex", "pop", " gws") ) {
		if ( var %in% names )
			if ( ! is.integer(phdata[[var]]))
				stop("variable ", var, " must be integer !")
	}

}

#' Convert a gwaa dataset to a GWS dataset and perform some sanity checks
#'
#' GWASBinTests datasets are GenABEL datasets with some special variables in the phenotypic variable data,
#' that we call GWS datasets.
#' The main variables are:
#' \describe{
#' 	\item{ pop }{ the Population index - It is the phenotype, for example case or control.
#'   You can currently have up to 256 different populations. This has to be coded as integers and be
#' 	\eqn{\le 255}
#'
#' 	Some functions such as \code{\link{processGws}} use this variable as the default phenotype, and
#'  require the all the phenotypes have to be defined, i.e. no NA is allowed. If you have NAs in your data
#'  you may use the \code{assignPopNATo} parameter to assign those NAs either to an existing phenotype or to
#'  a new one. Otherwise you may subset your Gws dataset to only consider the individuas with defined phenotypes.
#'
#' }
#'
#' \item{ gws }{ the Genome Wide Study index - It is used to perform \emph{meta-analysis} of several studies,
#' 	meaning that the permutations of the labels will be performed intra studies, i.e. a permutation will swap
#' 	labels from samples inside a same study. See \code{\link{mergeGws}} for merging studies in a GWS dataset.
#' 	This index has to be an integer.
#' 	}
#'
#' }
#'
#' @param gwaa The \code{\link[GenABEL]{gwaa.data-class}} object to convert
#'
#' @param assignPopNATo if defined, the NAs in the pop variable (or from the \var{bt} GenABEL variable if pop is not present)
#' 		will be assigned to this value.
#'
#' @return A \code{\link[GenABEL]{gwaa.data-class}} object converted. In fact only the phenodata \strong{sex} and \strong{gws}
#' 		might be modified.
#'
#'
#' @examples
#' 	data(srdta)
#'  # make a GWS dataset by setting the pop to 0
#' 	gws <- asGws( srdta, assignPopNATo=0)
#' @export
#'
#'
asGws <- function(gwaa, assignPopNATo=NULL) {

	if ( is.null(gwaa@phdata$pop) ) {
		if ( ! is.null(gwaa@phdata$bt) )
			gwaa@phdata$pop <- as.integer(gwaa@phdata$bt)
		else
			stop("no phenotype pop found")
	}

	if ( ! is.null(assignPopNATo)  )
		gwaa@phdata$pop[is.na(gwaa@phdata$pop)] <- as.integer(assignPopNATo)

	if ( any( is.na(gwaa@phdata$pop)) )
		stop("no NA pop phenotype allowed")

	gwaa@phdata$sex <- as.integer(gwaa@phdata$sex)
	if ( ! is.null(gwaa@phdata$gws) )
		gwaa@phdata$gws <- as.integer(gwaa@phdata$gws)

	return(gwaa)
}


#' Read a genome wide association study dataset in PLINK transposed-ped format
#'
#' This function is a quite easy way to read and convert a PLINK dataset.
#' It still requires that you convert the PLINK dataset to transposed-ped format (.tped and .tfam).
#'
#' If you have the PLINK software installed, you just have to use the \code{--recode} and \code{--transpose options}.
#' For example, if you have a binary PLINK dataset named foo (foo.bed, foo.bim, foo.fam), you can
#' convert it into the bar transposed dataset:
#'
#' \code{plink --recode --transpose --bfile foo --out bar}
#'
#' The dataset should be stored in two files named \code{basename}\strong{.gag}
#'  and \code{basename}\strong{.gap} for the genotypes and the phenotypes.
#'
#' @param basename Full path to files without extensions.
#'   The function adds successively \code{.tped} and \code{.tfam} to the path in order to have the full paths to files.
#' @param verbose if TRUE, suppress the output of the GenABEL \code{\link[GenABEL]{load.gwaa.data}} function.
#'
#' @return A genome wide association dataset as a \code{\link[GenABEL]{gwaa.data-class}} object.
#'
#' @export
#'
readPlinkTransposedData<- function(basename, verbose = FALSE)
{
	GENABEL_GENOTYPES_SUFFIX <- '.gag'
	GENABEL_PHENOTYPES_SUFFIX <- '.gap'

	outfile <- tempfile("outfile_read_plink_tped")

	genofile <- paste(outfile,GENABEL_GENOTYPES_SUFFIX, sep="")
	phenofile  <- paste(outfile,GENABEL_PHENOTYPES_SUFFIX, sep="")
	tpedfile  <- paste(basename,".tped", sep="")
	tfamfile <- paste(basename,".tfam", sep="")

	load <- function() convert.snp.tped(tpedfile, tfamfile, genofile)

	if ( verbose )
		load()
	else
		output <- capture.output( gws <- load() )

	# == read and convert tfam to phenofile
	# 6 columns: Family ID , Individual  id, father, mother, sex(1=male, 2=female, other=unknown, disease

	tfam <- read.delim(tfamfile, header = FALSE, stringsAsFactors = FALSE, sep = ' ')
	names(tfam) <- c("familyID", "id", "father", "mother", "sex", "pop")
	pheno <- tfam[c("id", "sex", "pop")]
	# convert sex to (1==male, females==0)
	pheno$sex <- pheno$sex %% 2

	write.table(pheno, file = phenofile, sep = ' ', row.names = FALSE)

	# now read it
	gws <- readGws(outfile)

	# now delete the files
	unlink(c(genofile, phenofile))

	.checkPhenoData(gws@phdata)

	return(gws)
}

#' Read a genome wide association study dataset (genotypes + phenotypes)
#'
#' Reads a genome wide association study dataset in the internal GenABEL format.
#' It is just a proxy for the \code{\link[GenABEL]{load.gwaa.data}} function of \code{\link[GenABEL]{GenABEL}}.
#'
#' The dataset should be stored in two files named \code{basename}\strong{.gag}
#'  and \code{basename}\strong{.gap}
#' for the genotypes and the phenotypes.
#'
#' @param basename Full path to files without extensions.
#'   The function adds successively \code{.gag} and \code{.gap} to the path in order to have the full paths to files.
#' @param verbose if TRUE, suppress the output of the GenABEL \code{\link[GenABEL]{load.gwaa.data}} function.
#'
#' @return A genome wide association dataset as a \code{\link[GenABEL]{gwaa.data-class}} object.
#'
#' @export
#'

readGws <- function(basename, verbose = FALSE)
{
	GENABEL_GENOTYPES_SUFFIX <- '.gag'
	GENABEL_PHENOTYPES_SUFFIX <- '.gap'

	load <- function() load.gwaa.data(phe = paste(basename,GENABEL_PHENOTYPES_SUFFIX, sep=""),
				gen = paste(basename,GENABEL_GENOTYPES_SUFFIX, sep=""))

	if ( verbose )
		gws <- load()
	else
		output <- capture.output( gws <- load() )

	.checkPhenoData(gws@phdata)

	return(gws)
}


#' Read a file containing defintions of genomic bins
#'
#' 	read a bins text file (TAB separated) into a dataframe,
#'   and sort the dataframe by genomic position.
#'
#' Bins are just genomic intervals, i.e. fully qualified by a chromosome and an
#' interval [start-end'].
#'
#' Bins can of course be restricted to a single location by defining start=end=position.
#' You could do this to ensure that some of the markers are alone in a bin.
#'
#'
#' \strong{FORMAT:}
#'   This file must be TAB separated, with a header line,
#' 		and contain the following fields:
#' 		\itemize{
#' 			\item{\code{chr_name}}{ name of chromosome (1-22, X, Y, ...);}
#' 			\item{\code{start}}{ start position of bin}
#' 			\item{\code{end}}{ end position of bin.}
#' 		}
#'
#' @param filename The path to the bins file.
#'
#' @return 	The sorted bins dataframe, with an additional column
#' 		\strong{bin_index}, useful to identify the bin
#'
#' @export
readBins <- function(filename)
{
	binDataFrame <- read.table(filename,header=TRUE,sep="\t")

	if (!("chr_name" %in% colnames(binDataFrame)))
		stop("ERROR: chr_name missing in ",filename)
	if (!("start" %in% colnames(binDataFrame)))
		stop("ERROR: start missing in ",filename)
	if (!("end" %in% colnames(binDataFrame)))
		stop("ERROR: end missing in ",filename)


#	# add the chromosome index column
	chrs <- convertChromosomeNames(binDataFrame$chr_name)

	# add a bin_indice column to the data frame
	nb_bins <- dim(binDataFrame)[1]
	binDataFrame$bin_index <- 1:nb_bins

	# make sure the file is sorted
	b.o <- order(chrs, binDataFrame$start)

	return(binDataFrame[b.o,]) ;
}

#' convert GenABEL chromosome names to integer
#'
#' Take a vector of chromosome names such as those use in a GenABEL dataset (\code{\link{chromosome}}
#' and convert them to integer.
#'
#' The non trivial conversions are: X=23, Y=24, XY=25, MT=26, NOTON(non localized)=0, MULTI(multi-localized)=-1
#' , UN(Unmapped)=-2, PAR(Pseudo Autosomal Region)=25, BAD(other error)=-100
#'
#' \strong{N.B:} In case of unknown name, the function will stop with an error.
#'
#' @param chr_names a vector of chromosome names or integers
#'
#' @return a vector of chromosome integer codes
#'
#' @examples
#' chrs <- convertChromosomeNames(c("1", "6", "x", "Y", "xY"))
#'
#' @export
convertChromosomeNames <-function(chr_names) {

	CHROMOSOME_NAMES_TO_INT <- list(X=23, Y=24,
			XY=25, # plink : Pseudo Autosomal Region
			MT=26, # plink: mitochondrial
			NOTON=0, # from dbsnp, mean no position found
			MULTI=-1,# from dbsnp, mean multi-localized
			UN=-2,# from dbsnp, mean localized but not on an assembly ?
			PAR=25,# dbsnp : Pseudo Autosomal Region
			BAD=-100) # Other reason

	chrs <- sapply(chr_names, function(chr_name) {

		if (is.integer(chr_name)) return(chr_name);

		i <- suppressWarnings( as.integer(as.character(chr_name)) )
		if ( ! is.na(i) )
			return(i);
		i <- CHROMOSOME_NAMES_TO_INT[[toupper(chr_name)]]
		if (is.null(i))
			stop(paste("unknown chromosome name:", chr_name))
		return(i)
	})
	return(as.integer(chrs))
}