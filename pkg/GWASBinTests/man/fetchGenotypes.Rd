\name{fetchGenotypes}
\alias{fetchGenotypes}
\title{fetch the genotypes from a GenABEL dataset in our format...}
\usage{fetchGenotypes(gws, index)
}
\description{fetch the genotypes from a \strong{GenABEL} dataset in our format}
\value{the genotypes, an integer matrix with values meaning: \describe{
\item{homozygous1 (aa) }{0}
\item{heterozygous (Aa) }{1}
\item{homozygous2 (AA) }{2}
\item{missing (\eqn{\emptyset}) }{-1}
}}
\arguments{\item{gws}{the \strong{GenABEL} dataset}
\item{index}{the index of the SNP. If NULL (the default), will fetch genotypes for all the SNPs}
}
\examples{data(srdta)
# one index
g <- fetchGenotypes(srdta,1)

# index range
gs <- fetchGenotypes(srdta, 1:10)

# names
gs <- fetchGenotypes(srdta, c("rs150", "rs179"))

# all
gs <- fetchGenotypes(srdta)}

