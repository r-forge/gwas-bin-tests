\name{processFiles}
\alias{processFiles}
\title{run the c++ analysis engine on data files...}
\usage{processFiles(basename="", geno_file="", pheno_file="", bins_file="",
    params=parameters(), covariables)
}
\description{run the c++ analysis engine on data files}
\details{This function does not load the data into R, it just passes the parameters to the C++ code.
It means that you do not need to load any data in R, everything happens in the NNBC C++ engine.}
\value{a data frame with the following columns:
\item{bin}{the bin index}
\item{chr}{the chr code (as integer, e.g. 23 for chromosome X )}
\item{start}{the start position of the bin}
\item{end}{the end position of the bin}
\item{nb_snps}{the number of SNPs in this bin}
\item{nb_permutations}{the actual number of permutations used to compute the pvalues of this bin}
\item{score_\strong{type}}{ the score for the given \strong{type}. For likelihoods it is not the ratio but the
likelihood under the alternative hypothesis }
\item{pv_\strong{type}}{ the pvalue for the given \strong{type} }
\item{fdr_\strong{type}}{ the \acronym{FDR} for the given \strong{type} }}
\seealso{\code{\link{processTable}}, \code{\link{processGws}}}
\arguments{\item{basename}{the basename of the files.
".gag", ".gap" and ".bins" suffixes will respectively be appended to the basename
to form the data file names. Purely for convenience.}
\item{geno_file}{the \link{GenABEL} genotypes file.
See  \code{\link[GenABEL]{load.gwaa.data}} for a description of the format.}
\item{pheno_file}{the \link{GenABEL} phenotypes data file, also refer to \code{\link[GenABEL]{load.gwaa.data}}.}
\item{bins_file}{the bins definition file. See \code{\link{readBins}} for more information about the format.}
\item{params}{the NNBC parameters, see \code{\link{parameters}}}
\item{covariables}{an optional dataframe of covariables. If no covariables are given, and
that the dataset stored in the files contain a phenotypic column "gws" it will be used by default as covariable}
}
\examples{\dontrun{
data_path <- system.file("extdata", package = "GWASBinTests")
basename <- paste(data_path, "/ms1", sep="")
res <- processFiles(basename)
}}

