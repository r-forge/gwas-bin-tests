\name{mergeGws}
\alias{mergeGws}
\title{Merge several genome wide association studies into a single dataset...}
\usage{mergeGws(gwsA, gwsB, intersected_snps_only=TRUE)
}
\description{Merge several genome wide association studies into a single dataset}
\details{Merge two association studies datasets as \code{\link[GenABEL]{gwaa.data-class}} objects into a single one.
The 'gws' column in phenodata keeps track of the sample study of origin.

\strong{N.B} the sample ids have to be distinct !!!

\strong{N.B} When using intersected_snps_only=FALSE, genotypes for samples for SNPs not in dataset
will be set to NA.



When doing a GWASBinTests analysis, the method takes into account the study of origin via the \strong{gws}
variable by permuting the phenotypes of the labels intra-study.

\strong{N.B}: the function is named \code{mergeGws} and not \code{merge.gws} to avoid \code{R CMD check}
warnings about \emph{S3 generic/method consistency} and the generic \code{\link{merge}} S3 function.}
\value{A genome wide association dataset as a \code{\link[GenABEL]{gwaa.data-class}} object.
Only the following fields are merged:

\item{pop}{ Case/Control status;}
\item{sex}{ Sex;}
\item{gws}{ Genome wide study identifier.}

If the \code{gws} field is absent in one of the studies,
it is added with a value different from the ones found in the other study.}
\arguments{\item{gwsA}{First genome wide association study as a \code{\link[GenABEL]{gwaa.data-class}} object.}
\item{gwsB}{Second genome wide association study as a \code{\link[GenABEL]{gwaa.data-class}} object.}
\item{intersected_snps_only}{merge only SNPs shared by gwsA and gwsB}
}

