\name{tryInterleaving}
\alias{tryInterleaving}
\alias{tryInterleavingSeq}
\alias{tryInterleavingThr}
\title{
  Test function for the package rngOpenMP.
}
\description{
  Test function for the package rngOpenMP in a sequential mode.
}
\usage{
  tryInterleavingSeq(ind)
  tryInterleavingThr(ind, nbThreads)
}
\arguments{
  \item{ind}{ Sequence of indices of the streams, from which the numbers
   should be generated. The indices are expected from \code{1:n} and
   are internally converted to C indexing \code{0:(n-1)}.}
  \item{nbThreads}{ The number of threads set using the function
        \code{omp_set_num_threads()} in the compiled code.}
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
  tryInterleavingSeq(ind)
  set.seed(12345)
  tryInterleavingThr(ind)
}

