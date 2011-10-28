#' Bin-based analysis of genome-wide association studies
#'
#' Given an a priori partitioning of the genome into regions termed \strong{bins},
#' GWASBinTests can compute p-values (and FDRs) of several association tests. Some are likelihood-based and
#' use genotyping error models that take into account genotyping errors and missing data,
#' some can take into account the correlation pattern of the neighbor markers, and some
#' are based on sum-scores.
#'
#' To be more precise, the available tests are:
#' \describe{
#' 		\item{ univariate }{ Univariate Genotypic Likelihood, corresponding to the \strong{L3} score
#' 			in the article. It is computed on genotype contigency tables, and it does not
#' 			take into account the correlation between markers of the same \strong{bin}.
#'  	}
#'
#' 		\item{ divariate }{ Divariate Genotypic Likelihood, corresponding to the \strong{L2} score
#' 			in the article.  It is computed on genotype contigency tables, and use the
#' 			correlation between consecutive pairs of markers.  }
#'
#' 		\item{ allelic }{ Allelic SumScore test. It is based on the sum of the allelic
#' 			(i.e. computed on the allelic contingency tables) pearson scores
#' 			of all the markers of a bin.
#' 			It should be very close to the set-based test implemented in \strong{PLINK} when
#' 			using the parameters \code{--set-max 99999 --set-p 1 --set-r2 1} }
#'
#' 		\item{ genotypic }{ Genotypic SumScore test. Like the allelic one but computed on genotypic tables. }
#' }
#'
#' \subsection{usage}{
#' 		GWASBinTests data are based on \code{\link{GenABEL}} datasets, and are called \strong{gws}.
#' 		The genotypes are stored using the gwaa-dataclass
#' 		and saved in \strong{raw} files, and the phenotype data in a data frame stored in a tabulated file.
#' 		GWASBinTests C++ engine can either be run directly on those data files, or on R data.
#'
#' 		The first step is to get and prepare some data. You can use samples from GenABEL (\code{\link{srdta}})
#' 		or directly from GWASBinTests (\code{\link{ms1}}).
#' 		You can convert data sets from PLINK format using \code{\link{readPlinkTransposedData}}, or convert genotypes using the
#' 		functions provided by \code{\link{GenABEL}}.
#' 		If you use a dataset from \code{\link{GenABEL}}, you have to adapt it, basically to set some phenotypic variables
#' 		that are used by GWASBinTests (see \code{\link{asGws}}).
#'
#' 		Then you need a bins description file (see \code{\link{readBins}}). There are is one available from the data/ directory
#' 		of GWASBinTests: \code{system.file("data/ms1.bins", package = "GWASBinTests")}
#'
#'		Now you can run the GWASBinTests analysis, using \code{\link{processFiles}} or \code{\link{processGws}}. Please take a look
#' 		at the parameters (see \code{\link{parameters}}), that can be fined tuned for performance and accuracy.
#'
#'		In the end you obtain a dataframe with all the results.
#'
#' 		You may also \emph{play} with hand-made tables to test the different p-values using \code{\link{processTable}}
#'
#' }
#'
#' \subsection{Implementation}{
#'
#' 		GWASBinTests is developped for reliability and efficiency.
#' 		To achieve these goals we developed two implementations of the tests.
#' 	    The first is optimized for simplicity, readability and conciseness and coded purely in R.
#'      It serves as a reference implementation, exhaustively tested with the included unit tests.
#' 		The second implementation is optimized for performance and thus coded in C++ and integrated in R via
#' 		\strong{Rcpp}.
#'
#' 		The pure-R implementation is too slow to be used on a big dataset, but can be used to validate the c++
#' 		optimized implementation.
#' 		There is also a standalone C++ executable not distributed for now but that available on request.
#' 	 	The integration with R has been greatly facilitated thanks to \code{\link{Rcpp}}.
#' 		We make use of the \strong{BOOST} C++ libraries, especially of the Binomial Distribution implementation
#' 		of confidence intervals.
#'
#'      We put a lot of efforts into the computational speed of the analysis,: \enumerate{
#' 			\item by carefully optimizing the bottlenecks of the code, especially the Random Number Generator and
#' 			the computation of the contingency tables
#'
#' 			\item by using multiple \strong{threads}. We used \strong{OpenMP} to implement the multihreading, because
#' 			it is easy and allow to write code that is both sequential and multithreaded. OpenMP should be available
#' 			with GCC \eqn{\ge 4.2}. The R package should detect automatically if OpenMP is available and compile the code
#' 			accordingly. See the \code{threads} parameter in \code{\link{parameters}} to control the number of threads that are used.
#'
#' 			\item by implementing \strong{heuristics} on the number of permutations to use to avoid useless rounds of computation.
#' 			See \code{\link{parameters}} for a description of these heuristics.
#'     }
#'
#' 		One of the bottlenecks was the Random Number Generator, both for speed and for the cycle length because
#' 		we may use a lot of permutations. Moreover we wanted to have reproducible results, even when multithreading.
#' 		All of these objectives have been met by using a modified version of the \strong{CRandomSFMT0} generator of the \strong{Randomc}
#' 		c++ library, which is \dQuote{an improvement of the Mersenne Twister with better randomness and
#' 		higher speed, designed specifically for processors with Single-Instruction-Multiple-Data
#' 		(SIMD) capabilities, such as the SSE2 and later instruction set}.
#' 		We further modified it by inlining some functions to achieve yet more speed.
#'		The cycle length is  \eqn{ \ge 2^{11213}-1} that allows a lot of independent permutations.
#' 		An unusual feature is that we get the very exact same results (at a given
#' 		fixed seed) in sequential or multithreaded modes, and thus obtain reproducible results, using multiple
#' 		seeds thanks to \emph{Randomc}.
#'
#'
#' 	   The result is that now we can compute in minutes what took days before.
#'
#' }
#'
#'
#'
#' @name GWASBinTests-package
#' @aliases GWASBinTests, NNBC
#' @docType package
#' @title p-values and FDRs for genomic regions in genome-wide association studies
#' @references
#' N.B: An application note for this package is in preparation.
#'
#' Articles presenting the method:
#' \describe{
#' \item{Long version:}{ \url{http://videolectures.net/site/normal_dl/tag=8590/article10.pdf} }
#' \item{Short version:}{ \url{http://www.biomedcentral.com/1753-6561/2/S4/S6} }
#' }
#' Additional resources: \url{http://videolectures.net/msht07_omont_gbba/}
#'
#' @author Nicolas Omont, Jerome Wojcik, Karl Forner
#' @seealso
#' \describe{
#'
#' \item{ PLINK }{ \url{http://pngu.mgh.harvard.edu/purcell/plink/}, \cite{Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
#' Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007)
#'         PLINK: a toolset for whole-genome association and population-based
#' 	linkage analysis. American Journal of Human Genetics, 81.}
#' Dirk Eddelbuettel, Romain Francois, with contributions by Simon}
#'
#' \item{ BOOST }{ free peer-reviewed portable C++ source libraries: \url{http://www.boost.org/} }
#'
#' \item{Randoma}{ package of random number generators: \url{http://www.agner.org/random/ran-instructions.pdf} }
#'
#' \item{OpenMP}{API specification for parallel programming: \url{http://openmp.org} }
#'
#'
#' }
#'
#' @keywords package
#'
NULL