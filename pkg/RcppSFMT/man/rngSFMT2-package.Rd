\name{rngSFMT2-package}
\alias{rngSFMT2-package}
\alias{rngSFMT2}
\docType{package}
\title{
  A R package providing C/C++/R interface to the SFMT Pseudo Randomn Generator based on the original implementation
}
\description{
 A R package providing C/C++/R interface to the SFMT Pseudo Randomn Generator based on the original implementation
}

\details{
The SFMT is a very fast RNG, so suitable for High Performance Computing tasks, which by nature
must be implemented in C/C++. This package provides an implementation of this RNG based on the original implementation
and provides R/C/C++ interfaces to help R developers using this RNG for their purposes.
The implementation is suitable for use in multi-threaded code, e.g. with OpenMP.

}


\seealso{
\describe{

\item{ SFMT }{the SFMT homepage, home of the original implementation: \url{http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/} }

\item{ Random Numbers in R }{Comments on random numbers in R \url{http://www2.cs.cas.cz/~savicky/randomNumbers/}}
}

}

\author{
  Petr Savicky
  Karl Forner

  Maintainer: Karl Forner <karl.forner@gmail.com>
}

\keyword{ package }

