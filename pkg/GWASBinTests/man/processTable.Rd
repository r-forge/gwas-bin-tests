\name{processTable}
\alias{processTable}
\title{run the c++ analysis engine on a univariate or divariate table...}
\usage{processTable(table, params=parameters())
}
\description{run the c++ analysis engine on a univariate or divariate table}
\details{Useful to study the different types of pvalues on some example tables.}
\value{a list :
\item{nb_permutations}{the actual number of permutations used to compute the pvalues of the table}
\item{pvalues}{a named vector of the pvalues computed for this table}}
\seealso{\code{\link{processFiles}}, \code{\link{processGws}}}
\arguments{\item{table}{a matrix of dimension \strong{nb_phenotype}*4 or \strong{nb_phenotype}*4*4
where table[i,j] (resp table[i,j,k] is the
number of samples for phenotype i and genotype j (resp phenotype i and genotypes (j,k).
The genotypes values are:
\enumerate{
\item{homozygous1}
\item{heterozygous}
\item{homozygous2}
\item{missing}
}}
\item{params}{the GWASBinTests parameters, see \code{\link{parameters}}}
}
\examples{table2x4 <- matrix(c(100, 0, 20, 1,  80, 10, 10, 0), nrow=2,ncol=4, byrow=TRUE)
res <- processTable(table2x4)

table2x4x4 <- c(277,250,5,2,0,0,0,1,
107,106,3,3,0,0,0,0,
8, 12, 1, 1,0,0,0,0,
1,1,0,0,0,0,0,0)
dim(table2x4x4) <- c(2,4,4)
res <- processTable(table2x4x4)}

