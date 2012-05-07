\name{trySetRngParams}
\alias{trySetRngParams}
\title{
Setting the parameters of generating streams of random numbers. Currently has
no effect.
}
\description{
Setting the parameters of generating streams of random numbers. The R level
interface for this function is meant to be used together with \code{getStreams}
for debugging.
}
\usage{
trySetRngParams(initType)
}
\arguments{
  \item{initType}{ The initialization type. Currently only 0 is supported, which 
   is the default.}
}
\details{
In future, the function may become used and may get further arguments.
}
\value{
  \code{NULL} invisibly.
}
\author{
Petr Savicky
}
\examples{
trySetRngParams(0)
}

