\name{fetchGenotypesAsList}
\alias{fetchGenotypesAsList}
\title{fetch the genotypes from a GenABEL dataset in our format as list of vectors...}
\usage{fetchGenotypesAsList(gws, index)
}
\description{fetch the genotypes from a \strong{GenABEL} dataset in our format as list of vectors}
\details{See \code{\link{fetchGenotypes}}}
\value{the genotypes, a list of genotypes integer vectors}
\arguments{\item{gws}{the \strong{GenABEL} dataset}
\item{index}{the index of the SNP(s). If NULL (the default), will fetch genotypes for all the SNPs}
}

