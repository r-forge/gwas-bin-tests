\name{generSFMT}
\alias{generIntSFMT}
\alias{generRealSFMT}
\alias{initSFMT}
\alias{deleteSFMT}
\title{Handling a SFMT object and generate random numbers from it.}
\description{
Create or delete a SFMT object in C++ and generate random numbers from it.
}
\usage{
initSFMT(seed=NULL, vseed=NULL)
generIntSFMT(n)
generRealSFMT(n)
deleteSFMT()
}
\arguments{
  \item{n}{ The required length of the output sequence.}
  \item{seed}{ Numeric value representing a single integer seed.}
  \item{vseed}{ A numeric vector represented an array to be used
     for initialization. Required, if \code{seed} is not specified.}
}
\details{
The parameters of SFMT are determined at compile time by 
"rngSFMT/src/Makevars" file.
}
\value{
\code{initSFMT()} and \code{deleteSFMT()} return \code{NULL}
invisibly. \code{generIntSFMT()} returns a numeric vector of
length \code{n} containing random 32-bit integer values.
Similarly, \code{generRealSFMT()} returns a numeric vector of
length \code{n} containing numbers from \code{(0, 1)}.
}
\examples{
initSFMT(vseed=1:3)
system.time(x <- generRealSFMT(1e7))
deleteSFMT()
}

