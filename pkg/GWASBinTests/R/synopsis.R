#' a sample session
#'
#' \subsection{data}{
#' 	Use \link{readGws}, \link{readPlinkTransposedData} and \link{readBins} to load/import your data.
#' 	See the \strong{GenABEL} documentation if you have data in other formats. Once you have a GenABEL dataset
#'  you can easily adapt it using \link{asGws}.
#'  Save your data with \link{saveGws}.
#' }
#'
#' A sample session to get you started. Be sure to understand the parameters, especially
#' those related to the heuristics, which can save you hours of computing time
#' (cf \strong{\link{parameters}})
#'
#'
#' @examples \dontrun{
#'
#' gws <- readGws("path/to/your_dataset")
#' bins <- readBins("path/to/your_bins_file.bins")
#'
#' params <- parameters(verbosity=1, nb_permutations=1000000, min_pvalue=0.01, max_relative_error=0.1)
#' results <- processGws(gws, bins, params=params, covariables=your_covariables)
#' }
#'
#' @name 0_Synopsis
#' @aliases synopsis
#' @title Synopsis
#' @keywords synopsis
NULL
