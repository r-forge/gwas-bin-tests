\name{testSums}
\alias{testSequential}
\alias{testThreaded}
\title{
Test function for the package rngOpenMP.
}
\description{
  Test and example functions for the package rngOpenMP using sequential
  and multithreaded mode under OpenMP.
}
\usage{
  testSequential(n, len)
  testThreaded(n, len)
}
\arguments{
  \item{n}{ Number of parallel streams required.}
  \item{len}{ A numeric vector of length \code{n} specifying the length
        of the sequence of numbers generated from each stream.}
}
\details{
  Generates \code{len[i]} numbers from stream \code{i} and computes their
  sum in the C code. The index \code{i} runs through \code{1:n}. The
  functions \code{testSequential} and \code{testThreaded} differ
  only in that the former runs sequentially and the latter uses OpenMP
  to parallelize the for loop over \code{i}.
}
\value{
  Named vector of the obtained sums of length \code{n}, where the names
  are the indices of the streams.
}
\author{
  Petr Savicky.
}
\examples{
  set.seed(12345)
  testSequential(5, c(16, 8, 4, 2, 1))
  set.seed(12345)
  testThreaded(5, c(16, 8, 4, 2, 1))
}

