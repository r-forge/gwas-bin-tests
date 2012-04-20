\name{testSFMT}
\alias{testSFMT11213}
\alias{testSFMT19937}
\alias{simulateSFMT11213}
\alias{simulateSFMT19937}
\title{R level simulation of SFMT suitable for test purposes.}
\description{
Functions \code{simulateSFMT11213()} and \code{simulateSFMT19937()}
provide a slow simulation of SFMT for the corresponding parameters using
an R script, which uses 0,1 vectors instead of integers represented in
binary.
}
\usage{
simulateSFMT11213(state, n)
simulateSFMT19937(state, n)
testSFMT11213(x)
testSFMT19937(x)
}
\arguments{
  \item{state}{ Numeric vector containing 32-bit integers representing
  an initial state of SFMT for the simulation.}
  \item{n}{ The length of the required output sequence obtained by
   simulation of SFMT by an R level script.}
  \item{x}{ Numeric vector containing 32-bit integers representing
  a sequence obtained from a C level SFMT.}
}
\details{
The functions \code{testSFMT11213(x)} and \code{testSFMT19937(x)} use an
initial part of \code{x} as the initial state of SFMT. The length of
this initial segment is 352 for SFMT11213 and 624 for SFMT19937. Then
the corresponding R level simulation function for this initial state
is called. The obtained output is compared to the rest of the sequence
\code{x} of length \code{length(x) - k}, where \code{k == 352} for
SFMT11213 and \code{k == 624} for SFMT19937.
}
\value{
Table of the frequencies of \code{TRUE} and \code{FALSE} obtained
by comparisons of individual components of the output sequence.
The total number of comparisons is \code{length(x) - k} as above.
}
\examples{
initSFMT(vseed=1234)
x <- generIntSFMT(1500)
deleteSFMT()
# The correct test depends on setting "MEXP" in the "rngSFMT/src/Makevars",
# when the package is compiled.
testSFMT19937(x)
testSFMT11213(x)
}

