\name{data#plink_small}
\alias{data#plink_small}
\alias{data#plink_small}
\title{A sample data set}
\description{This is a subset of an example dataset available from the Plink website}
\details{It has been downloaded and converted from the the files wgas1.ped and wgas1.map from the example.zip file.

The command used to generate it is:

\code{plink --chr 21 --from-kb 5000 --to-kb 20000  --recode --transpose --file wgas1 --out plink_small}


It consists of:
\itemize{
\item{ plink_small.tped }{ the transposed genotypes}
\item{ plink.tfam }{ the phenotypes}
}}
\alias{plink_small}
\docType{data}
\references{url{http://pngu.mgh.harvard.edu/~purcell/plink/res.shtml#teach}}
\keyword{data}

