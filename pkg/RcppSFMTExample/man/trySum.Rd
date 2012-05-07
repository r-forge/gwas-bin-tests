\name{trySum}
\alias{trySum}
\alias{trySumSeq}
\alias{trySumThr}
\title{
Test function for the package rngOpenMP.
}
\description{
  Test and example functions for the package rngOpenMP using sequential
  and multithreaded mode under OpenMP.
}
\usage{
  trySumSeq(len)
  trySumThr(len, nbThreads)
}
\arguments{
  \item{len}{ A numeric vector of length \code{n} specifying the length
        of the sequence of numbers generated from each stream.}
  \item{nbThreads}{ The number of threads set using the function
        \code{omp_set_num_threads()} in the compiled code.}
}
\details{
  Generates \code{len[i]} numbers from stream \code{i} and computes their
  sum in the C code. The index \code{i} runs through \code{1:n}. The
  functions \code{trySumSeq} and \code{trySumThr} differ
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
  trySumSeq(c(16, 8, 4, 2, 1))
  set.seed(12345)
  trySumThr(c(16, 8, 4, 2, 1))
}

