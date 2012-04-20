\name{formatOutput}
\alias{formatOutput}
\title{Prints 1000 random integers for test purposes.}
\description{
Prints 1000 random integers in the format suitable for comparison to
the precomputed values contained in the distribution files of SFMT.
}
\usage{
formatOutput(useArray)
}
\arguments{
  \item{useArray}{ Logical, whether the test numbers for the single integer
  initialization or for the initialization by array should be printed.}
}
\details{
The seeds used in both variants are the same as in "test.c" file
in SFMT distribution.
}
\value{
\code{NULL} invisibly.
}
\examples{
# Compare the created file "out.txt" with "SFMT.11213.out.txt" or "SFMT.19937.out.txt"
# from the distribution of SFMT. The choice of "MEXP" depends on its setting in
# "rngSFMT/src/Makevars" file, when the package is compiled.
\dontrun{
sink("out.txt")
formatOutput(useArray=FALSE)
formatOutput(useArray=TRUE)
sink()
}
}

