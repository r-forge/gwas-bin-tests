\name{convertChromosomeNames}
\alias{convertChromosomeNames}
\title{convert GenABEL chromosome names to integer...}
\usage{convertChromosomeNames(chr_names)
}
\description{convert GenABEL chromosome names to integer}
\details{Take a vector of chromosome names such as those use in a GenABEL dataset (\code{\link{chromosome}}
and convert them to integer.

The non trivial conversions are: X=23, Y=24, XY=25, MT=26, NOTON(non localized)=0, MULTI(multi-localized)=-1
, UN(Unmapped)=-2, PAR(Pseudo Autosomal Region)=25, BAD(other error)=-100

\strong{N.B:} In case of unknown name, the function will stop with an error.}
\value{a vector of chromosome integer codes}
\arguments{\item{chr_names}{a vector of chromosome names or integers}
}
\examples{chrs <- convertChromosomeNames(c("1", "6", "x", "Y", "xY"))}

