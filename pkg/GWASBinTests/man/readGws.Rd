\name{readGws}
\alias{readGws}
\title{Read a genome wide association study dataset (genotypes + phenotypes)...}
\usage{readGws(basename, verbose=FALSE)
}
\description{Read a genome wide association study dataset (genotypes + phenotypes)}
\details{Reads a genome wide association study dataset in the internal GenABEL format.
It is just a proxy for the \code{\link[GenABEL]{load.gwaa.data}} function of \code{\link[GenABEL]{GenABEL}}.

The dataset should be stored in two files named \code{basename}\strong{.gag}
and \code{basename}\strong{.gap}
for the genotypes and the phenotypes.}
\value{A genome wide association dataset as a \code{\link[GenABEL]{gwaa.data-class}} object.}
\arguments{\item{basename}{Full path to files without extensions.
The function adds successively \code{.gag} and \code{.gap} to the path in order to have the full paths to files.}
\item{verbose}{if TRUE, suppress the output of the GenABEL \code{\link[GenABEL]{load.gwaa.data}} function.}
}

