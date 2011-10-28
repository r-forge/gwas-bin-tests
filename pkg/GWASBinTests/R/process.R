# functions related to the C++ nnbc engine

# internal function used to format the data sent from the C++ Rcpp layer
# into a data frame
.convert_nnbc_results_to_data_frame <- function(res) {
	# transform pvalues_by_bin and fdrs_by_bin
	# from list of numeric vectors to matrices

	if (length(res) == 0)
		return(data.frame())

	nb_rows <- length(res$chr)
	types <- res$types
	nb_pvalues <- length(types)
	pvalues <- matrix(unlist(res$pvalues_by_bin), nrow = nb_rows, ncol = nb_pvalues, byrow = TRUE)
	scores <- matrix(unlist(res$scores_by_bin), nrow = nb_rows, ncol = nb_pvalues, byrow = TRUE)
	fdrs <- matrix(unlist(res$fdrs_by_bin), nrow = nb_rows, ncol = nb_pvalues, byrow = TRUE)

	colnames(pvalues) <- paste("pv_", types, sep="")
	colnames(scores) <- paste("score_", types, sep="")
	colnames(fdrs) <- paste("fdr_", types, sep="")

	df <- data.frame(res[c("bin", "chr", "start", "end", "nb_snps", "nb_permutations")], pvalues, scores, fdrs)

	return(df)
}

# internal function that converts a dataframe of covariable to a vector of group numbers
.compute_groups_from_covariables <- function(covariables) {
	if ( is.null(covariables) || nrow(covariables) == 0 )
		return(integer())

	nb <- nrow(covariables)
	uni <- unique(covariables)
	covariables$order <- 1:nb
	uni$key <- 1:nrow(uni)
	m <- merge(covariables, uni)
	groups <- m[order(m$order),]$key
	return(as.integer(groups))
}


#' run the c++ analysis engine on data files
#'
#' This function does not load the data into R, it just passes the parameters to the C++ code.
#' It means that you do not need to load any data in R, everything happens in the NNBC C++ engine.
#'
#' @param basename the basename of the files.
#' 	".gag", ".gap" and ".bins" suffixes will respectively be appended to the basename
#' 	 to form the data file names. Purely for convenience.
#' @param geno_file the \link{GenABEL} genotypes file.
#' 	See  \code{\link[GenABEL]{load.gwaa.data}} for a description of the format.
#' @param pheno_file the \link{GenABEL} phenotypes data file, also refer to \code{\link[GenABEL]{load.gwaa.data}}.
#' @param bins_file the bins definition file. See \code{\link{readBins}} for more information about the format.
#' @param params the NNBC parameters, see \code{\link{parameters}}
#' @param covariables an optional dataframe of covariables. If no covariables are given, and
#' 		that the dataset stored in the files contain a phenotypic column "gws" it will be used by default as covariable
#'
#' @return a data frame with the following columns:
#' \item{bin}{the bin index}
#' \item{chr}{the chr code (as integer, e.g. 23 for chromosome X )}
#' \item{start}{the start position of the bin}
#' \item{end}{the end position of the bin}
#' \item{nb_snps}{the number of SNPs in this bin}
#' \item{nb_permutations}{the actual number of permutations used to compute the pvalues of this bin}
#' \item{score_\strong{type}}{ the score for the given \strong{type}. For likelihoods it is not the ratio but the
#' likelihood under the alternative hypothesis }
#' \item{pv_\strong{type}}{ the pvalue for the given \strong{type} }
#' \item{fdr_\strong{type}}{ the \acronym{FDR} for the given \strong{type} }
#'
#' @seealso \code{\link{processTable}}, \code{\link{processGws}}
#'
#' @examples
#' \dontrun{
#' data_path <- system.file("extdata", package = "GWASBinTests")
#' basename <- paste(data_path, "/ms1", sep="")
#' res <- processFiles(basename)
#' }
#'
#' @export

processFiles <- function(basename="", geno_file = "", pheno_file = "", bins_file ="",
		params = parameters(), covariables=NULL ) {

	#geno_file <- pheno_file <- bins_file <- ""
	if ( ! nchar( geno_file) )
		geno_file <- paste(basename,".gag",sep="")
	if ( ! nchar( pheno_file) )
		pheno_file <- paste(basename,".gap",sep="")
	if ( ! nchar( bins_file) )
		bins_file <- paste(basename,".bins",sep="")

	if(  ! ( nchar( geno_file) &&  nchar( pheno_file) && nchar( bins_file)) )
		stop("One of the data files is not defined")

	files <- list(geno_file=geno_file, pheno_file=pheno_file, bins_file=bins_file)

	parameters <- .checkParameters(params)
	groups <- GWASBinTests:::.compute_groups_from_covariables(covariables)

	res <- .Call("RcppRunNNBCOnFiles", files, parameters, groups,
			PACKAGE="GWASBinTests")

	return( .convert_nnbc_results_to_data_frame(res))
}



#' run the c++ analysis engine on a univariate or divariate table
#'
#' Useful to study the different types of pvalues on some example tables.
#'
#' @param table a matrix of dimension \strong{nb_phenotype}*4 or \strong{nb_phenotype}*4*4
#'  where table[i,j] (resp table[i,j,k] is the
#' 	number of samples for phenotype i and genotype j (resp phenotype i and genotypes (j,k).
#' The genotypes values are:
#' \enumerate{
#' 		\item{homozygous1}
#' 		\item{heterozygous}
#' 		\item{homozygous2}
#' 		\item{missing}
#' }
#'
#' @param params the GWASBinTests parameters, see \code{\link{parameters}}
#' @return a list :
#' \item{nb_permutations}{the actual number of permutations used to compute the pvalues of the table}
#' \item{pvalues}{a named vector of the pvalues computed for this table}
#'
#' @seealso  \code{\link{processFiles}}, \code{\link{processGws}}
#'
#' @examples
#' 	table2x4 <- matrix(c(100, 0, 20, 1,  80, 10, 10, 0), nrow=2,ncol=4, byrow=TRUE)
#'  res <- processTable(table2x4)
#'
#' 	table2x4x4 <- c(277,250,5,2,0,0,0,1,
#'  	107,106,3,3,0,0,0,0,
#' 		8, 12, 1, 1,0,0,0,0,
#' 		1,1,0,0,0,0,0,0)
#'  dim(table2x4x4) <- c(2,4,4)
#'  res <- processTable(table2x4x4)
#'
#' @export
processTable <- function(table, params = parameters()) {

	dims <- dim(table)
	nb_dims <- length(dims)

	if ( nb_dims < 2 || nb_dims > 3)
		stop("argument table must be a nx4 or nx4x4 matrix!")
	nb_phenotypes <- dims[1]
	if ( nb_phenotypes < 2)
		stop("nb_phenotypes, the first dimension of the matrix must be >= 2")

	parameters <- .checkParameters(params)

	v <- as.integer(table)
	res <- .Call("RcppRunNNBCOnTable", nb_phenotypes, v, parameters,
			PACKAGE="GWASBinTests")

	# res is a list with nb_permutations, pvalues, types

	nb_pvalues <- length(res$pvalues)
	names(res$pvalues) <- res$types

	return(res)
}

#' run the c++ analysis engine on a GenABEL gwaa dataset
#'
#'
#' The GenABEL dataset needs to be \emph{GWSified} in order to be processed,
#' see \code{\link{asGws}}
#'
#' @param gws a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}
#' @param bins a data frame, as returned by \code{\link{readBins}}. If no bins are given
#' 		the function will use a bin for each SNP to simulate a marker by marker scan
#' 		(see \code{\link{makeBinsForSnps}}
#' @param params the GWASBinTests parameters, see \code{\link{parameters}}
#' @param phenotypes an optional vector of phenotypes. It will be converted to integer and must not
#' 	contain a value higher than 255
#' @param covariables an optional dataframe of covariables. If no covariables are given, and
#' 		that the dataset stored in the files contain a phenotypic column \strong{gws} it will be used by default as covariable
#'
#' @return a data frame (cf \code{\link{processFiles}})
#'
#' @examples
#'	data("ms1")
#'  res <- processGws(ms1_gws, ms1_bins,
#'    parameters(nb_permutations=1000000, max_relative_error=0.1)
#'  )
#'
#' @export
processGws <- function(gws, bins=NULL, params = parameters(), phenotypes=phdata(gws)$pop, covariables=NULL) {
	gt <- gtdata(gws)
	# ensure chromosome are encoded as strings
	if ( is.null(bins ) )
		bins <- makeBinsForSnps(gws)
	bins[[1]] <- as.character(bins[[1]])

	parameters <- .checkParameters(params)

	phdata <- phdata(gws)[c("id", "sex", "pop")]
	phdata$sex <- as.integer(phdata$sex)


	groups <- GWASBinTests:::.compute_groups_from_covariables(covariables)
	if ( length(groups) == 0 && "gws" %in% names(gws@phdata) )
		groups <- as.integer(gws@phdata$gws)

	# check the phenotypes: no NA allowed
	phenotypes <- as.integer(phenotypes)
	if ( any(is.na(phenotypes)) )
		stop("no NA allowed in phenotypes")
	# check the max
	if ( max(phenotypes) > 255 )
		stop("Error, phenotype must be <= 255")
	# check the length
	if ( length(phenotypes) != nids(gws))
		stop("Bad number of phenotypes")

	# send it via the pop phdata column
	phdata$pop <- phenotypes

	res <- .Call("RcppRunNNBConRdata"
			,idnames(gt)
			,map(gt)
			,gt@coding # we prefer as integer than strings
			,snpnames(gt)
			,chromosome(gt)
			,gt@strand # we prefer as integer than strings
			,gt@gtps@.Data
			,phdata
			,bins,
			parameters,
			groups,
			PACKAGE="GWASBinTests"
	)
	return( GWASBinTests:::.convert_nnbc_results_to_data_frame(res))
}
