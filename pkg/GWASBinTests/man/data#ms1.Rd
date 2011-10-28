\name{data#ms1}
\alias{data#ms1}
\alias{data#ms1}
\title{A sample data set}
\description{This is the first 1000 SNPs of an internal and private dataset.}
\details{The sample names have been anonymized.


The Bins are defined on DNA from protein genes as defined in the version 35:35 of EnsEMBL Birney et al.
(2006) of the human DNA sequence. The basic region
of a gene lies from the beginning of its first exon
to the end of its last exon. Overlapping genes are clustered
in the same bin. If two consecutive genes or clusters
of overlapping genes are separated by less than
200 kbp, the bin limit is fixed in the middle of the interval.
Otherwise, the limit of the upstream bin is set
50 kbp downstream its last exon, the limit of the downstream
bin is set 50 kbp upstream its first exon, and a
special bin corresponding to a desert is created in between
the two bins. With these rules, desert bins have
a minimum length of 100 kbp

It consists of:
\itemize{
\item{ms1.gag}{ the genotypes}
\item{ms1.gap}{ the phenotypes}
\item{ms1.bins}{ the bins definitions}
}

It can be loaded via the data() function:
\code{data("ms1")}, and provides:
\describe{
\item{ \code{ms1_gws} }{ a \emph{GWSified} \code{\link[GenABEL]{gwaa.data-class}} dataset }
\item{ \code{ms1_bins} }{ a bins description dataframe, as returned by \code{\link{readBins}} }

}}
\alias{ms1}
\docType{data}
\keyword{data}

