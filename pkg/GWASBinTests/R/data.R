# data.R DESCRIBE datasets in data/

# for documenting data sets, the only way I found is to place
# the description first, without any tag

#' This is the first 1000 SNPs of an internal and private dataset.
#'
#' The sample names have been anonymized.
#'
#'
#' The Bins are defined on DNA from protein genes as defined in the version 35:35 of EnsEMBL Birney et al.
#' (2006) of the human DNA sequence. The basic region
#' of a gene lies from the beginning of its first exon
#' to the end of its last exon. Overlapping genes are clustered
#' in the same bin. If two consecutive genes or clusters
#' of overlapping genes are separated by less than
#' 200 kbp, the bin limit is fixed in the middle of the interval.
#' Otherwise, the limit of the upstream bin is set
#' 50 kbp downstream its last exon, the limit of the downstream
#' bin is set 50 kbp upstream its first exon, and a
#' special bin corresponding to a desert is created in between
#' the two bins. With these rules, desert bins have
#' a minimum length of 100 kbp
#'
#' It consists of:
#' \itemize{
#' 		\item{ms1.gag}{ the genotypes}
#' 		\item{ms1.gap}{ the phenotypes}
#' 		\item{ms1.bins}{ the bins definitions}
#' }
#'
#' It can be loaded via the data() function:
#' \code{data("ms1")}, and provides:
#' \describe{
#'  \item{ \code{ms1_gws} }{ a \emph{GWSified} \code{\link[GenABEL]{gwaa.data-class}} dataset }
#'  \item{ \code{ms1_bins} }{ a bins description dataframe, as returned by \code{\link{readBins}} }
#'
#' }
#'
#' @name data#ms1
#' @aliases ms1
#'
#' @title A sample data set
#' @docType data
#' @keywords data
NULL

#' See \link{data#ms1}
#' @name data#ms1_gws
#' @aliases ms1_gws
#' @docType data
#' @keywords data
NULL

#' See \link{data#ms1}
#' @name data#ms1_bins
#' @aliases ms1_bins
#' @docType data
#' @keywords data
NULL

#' This is a subset of an example dataset available from the Plink website
#'
#' It has been downloaded and converted from the the files wgas1.ped and wgas1.map from the example.zip file.
#'
#' The command used to generate it is:
#'
#' \code{plink --chr 21 --from-kb 5000 --to-kb 20000  --recode --transpose --file wgas1 --out plink_small}
#'
#'
#' It consists of:
#' \itemize{
#' 		\item{ plink_small.tped }{ the transposed genotypes}
#' 		\item{ plink.tfam }{ the phenotypes}
#' }
#' @name data#plink_small
#' @aliases plink_small
#' @title A sample data set
#' @docType data
#' @references url{http://pngu.mgh.harvard.edu/~purcell/plink/res.shtml#teach}
#' @keywords data
NULL

