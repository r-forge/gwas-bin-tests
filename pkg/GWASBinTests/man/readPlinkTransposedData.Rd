\name{readPlinkTransposedData}
\alias{readPlinkTransposedData}
\title{Read a genome wide association study dataset in PLINK transposed-ped format...}
\usage{readPlinkTransposedData(basename, verbose=FALSE)
}
\description{Read a genome wide association study dataset in PLINK transposed-ped format}
\details{This function is a quite easy way to read and convert a PLINK dataset.
It still requires that you convert the PLINK dataset to transposed-ped format (.tped and .tfam).

If you have the PLINK software installed, you just have to use the \code{--recode} and \code{--transpose options}.
For example, if you have a binary PLINK dataset named foo (foo.bed, foo.bim, foo.fam), you can
convert it into the bar transposed dataset:

\code{plink --recode --transpose --bfile foo --out bar}

The dataset should be stored in two files named \code{basename}\strong{.gag}
and \code{basename}\strong{.gap} for the genotypes and the phenotypes.}
\value{A genome wide association dataset as a \code{\link[GenABEL]{gwaa.data-class}} object.}
\arguments{\item{basename}{Full path to files without extensions.
The function adds successively \code{.tped} and \code{.tfam} to the path in order to have the full paths to files.}
\item{verbose}{if TRUE, suppress the output of the GenABEL \code{\link[GenABEL]{load.gwaa.data}} function.}
}

