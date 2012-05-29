\name{SFMT}
\alias{SFMT}
\title{
the SFMT class is a reference class wrapping the c++ class implementing the SFMT Random Number Generator
}
\description{
This class allows to use instances of the c++ RNG. 
type: SFMT for a brief description of the class, its constructors and methods.

}



\author{
Karl Forner
}
\examples{
r <- new(SFMT, 1234)
x <- r$nextRandomInteger()
}

