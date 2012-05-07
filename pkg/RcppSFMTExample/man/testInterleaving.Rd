\name{testInterleaving}
\alias{testInterleaving}
\title{
  Test function for the package rngOpenMP.
}
\description{
  Test function for the package rngOpenMP in a sequential mode.
}
\usage{
  testInterleaving(n, ind)
}
\arguments{
  \item{n}{ Number of parallel streams required.}
  \item{ind}{ Sequence of indices of the streams, from which the numbers
   should be generated. The indices are expected from \code{1:n} and
   are internally converted to C indexing \code{0:(n-1)}.}
}
\details{
  The function calls a C-level function in rngVerify, which calls
  C-functions from rngOpenMP to generate random numbers from \code{n}
  streams in a sequential mode. The output is a single stream of random
  numbers, where \code{i}-th number is from stream \code{ind[i]}. The 
  function tests that generating a number from one stream does not
  affect the implementation of other streams in a sequential mode.
}
\value{
  Named vector of the generated numbers. The names specify the stream,
  to which each number belongs.
}
\author{
  Petr savicky
}
\examples{
  n <- 3
  ind <- c(1, 2, 1, 3, 2, 1)
  set.seed(12345)
  testInterleaving(n, ind)
  set.seed(12345)
  testInterleaving(n, sort(ind))
}

