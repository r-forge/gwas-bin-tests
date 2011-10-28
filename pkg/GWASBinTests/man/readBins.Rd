\name{readBins}
\alias{readBins}
\title{Read a file containing defintions of genomic bins...}
\usage{readBins(filename)
}
\description{Read a file containing defintions of genomic bins}
\details{read a bins text file (TAB separated) into a dataframe,
and sort the dataframe by genomic position.

Bins are just genomic intervals, i.e. fully qualified by a chromosome and an
interval [start-end'].

Bins can of course be restricted to a single location by defining start=end=position.
You could do this to ensure that some of the markers are alone in a bin.


\strong{FORMAT:}
This file must be TAB separated, with a header line,
and contain the following fields:
\itemize{
\item{\code{chr_name}}{ name of chromosome (1-22, X, Y, ...);}
\item{\code{start}}{ start position of bin}
\item{\code{end}}{ end position of bin.}
}}
\value{The sorted bins dataframe, with an additional column
\strong{bin_index}, useful to identify the bin}
\arguments{\item{filename}{The path to the bins file.}
}

