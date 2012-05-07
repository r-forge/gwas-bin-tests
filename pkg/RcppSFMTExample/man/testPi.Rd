\name{testPi}
\alias{testPi}
\title{
Test function for the package rngOpenMP.
}
\description{
  Computes pi by counting how many random points in a square
  belong to a circle using calls of functions of OpenMP package
  in multithreaded mode under OpenMP.
}
\usage{
  testPi(n, nbThreads=8)
}
\arguments{
  \item{n}{ Number of sample points.}
  \item{nbThreads}{ The number of threads set using the function
        \code{omp_set_num_threads()} in the compiled code.}
}
\value{
  A numeric estimate of pi.
}
\examples{
  set.seed(12345)
  testPi(5e6)
}

