\name{asGws}
\alias{asGws}
\title{Convert a gwaa dataset to a GWS dataset and perform some sanity checks...}
\usage{asGws(gwaa, assignPopNATo)
}
\description{Convert a gwaa dataset to a GWS dataset and perform some sanity checks}
\details{GWASBinTests datasets are GenABEL datasets with some special variables in the phenotypic variable data,
that we call GWS datasets.
The main variables are:
\describe{
\item{ pop }{ the Population index - It is the phenotype, for example case or control.
You can currently have up to 256 different populations. This has to be coded as integers and be
\eqn{\le 255}

Some functions such as \code{\link{processGws}} use this variable as the default phenotype, and
require the all the phenotypes have to be defined, i.e. no NA is allowed. If you have NAs in your data
you may use the \code{assignPopNATo} parameter to assign those NAs either to an existing phenotype or to
a new one. Otherwise you may subset your Gws dataset to only consider the individuas with defined phenotypes.

}

\item{ gws }{ the Genome Wide Study index - It is used to perform \emph{meta-analysis} of several studies,
meaning that the permutations of the labels will be performed intra studies, i.e. a permutation will swap
labels from samples inside a same study. See \code{\link{mergeGws}} for merging studies in a GWS dataset.
This index has to be an integer.
}

}}
\value{A \code{\link[GenABEL]{gwaa.data-class}} object converted. In fact only the phenodata \strong{sex} and \strong{gws}
might be modified.}
\arguments{\item{gwaa}{The \code{\link[GenABEL]{gwaa.data-class}} object to convert}
\item{assignPopNATo}{if defined, the NAs in the pop variable (or from the \var{bt} GenABEL variable if pop is not present)
will be assigned to this value.}
}
\examples{data(srdta)
# make a GWS dataset by setting the pop to 0
gws <- asGws( srdta, assignPopNATo=0)}

