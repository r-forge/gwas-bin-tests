\name{processGws}
\alias{processGws}
\title{run the c++ analysis engine on a GenABEL gwaa dataset...}
\usage{processGws(gws, bins, params=parameters(), phenotypes=phdata(gws)$pop, covariables)
}
\description{run the c++ analysis engine on a GenABEL gwaa dataset}
\details{The GenABEL dataset needs to be \emph{GWSified} in order to be processed,
see \code{\link{asGws}}}
\value{a data frame (cf \code{\link{processFiles}})}
\arguments{\item{gws}{a GWSified dataset of class \code{\link[GenABEL]{gwaa.data-class}}}
\item{bins}{a data frame, as returned by \code{\link{readBins}}. If no bins are given
the function will use a bin for each SNP to simulate a marker by marker scan
(see \code{\link{makeBinsForSnps}}}
\item{params}{the GWASBinTests parameters, see \code{\link{parameters}}}
\item{phenotypes}{an optional vector of phenotypes. It will be converted to integer and must not
contain a value higher than 255}
\item{covariables}{an optional dataframe of covariables. If no covariables are given, and
that the dataset stored in the files contain a phenotypic column \strong{gws} it will be used by default as covariable}
}
\examples{data("ms1")
res <- processGws(ms1_gws, ms1_bins,
parameters(nb_permutations=1000000, max_relative_error=0.1)
)}

