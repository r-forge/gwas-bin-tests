\name{simpleGenotypeErrorModel}
\alias{simpleGenotypeErrorModel}
\title{compute the simple genotype error model for a marker (SNP)...}
\usage{simpleGenotypeErrorModel(alpha, beta)
}
\description{compute the \emph{simple} genotype error model for a marker (SNP)}
\details{See \link{concept#genotype.error.model}.

This model is much simpler than the affymetrix one, in that
in does not assume that there is a bias on the heterozygous genotypes.

This will generate the following genotype error model:

\deqn{
P( \mathcal{O}^j \backslash  \mathcal{G}^j ) =
\left(
\begin{array}{r|r|r|r}
\mathcal{O}^j \backslash  \mathcal{G}^j  & aa & Aa & AA\\
\hline

aa & (1-\beta)*(1-\alpha) & (1-\beta)*\alpha/2 &  (1-\beta)*\alpha/2 \\
\hline
Aa & (1-\beta)*\alpha/2 & (1-\beta)*(1-\alpha) & (1-\beta)*\alpha/2 \\
\hline
AA & 1-\beta)*\alpha/2 & 1-\beta)*\alpha/2 & (1-\beta)*(1-\alpha) \\
\hline
\emptyset &  \beta & \beta & \beta \\
\end{array}
\right)
}

Some comments: \enumerate{
\item \eqn{ (1-\beta) } is the probability for an observed genotype not to be missing.
That is why this term is included in all non-missing observed genotype rows.
\item \eqn{ (1-\alpha) } is the probability of not doing an error
\item the probability of not doing any mistake (the diagonal) for a given genotype is \eqn{(1-\beta)*(1-\alpha)},
which simply means that this is not missing AND there is no error.
\item the error is equally split among the two other genotypes

}}
\value{return error model a double}
\arguments{\item{alpha}{the error_rate \eqn{ \alpha }. See \strong{\link{parameters}}}
\item{beta}{the missing value rate \eqn{ \beta }}
}
\examples{m <- simpleGenotypeErrorModel(0.05, 0.01)

# error model for first SNP of srdta
# by computing missing rate on genotypes
data(srdta)
beta <- missingRateFromGenotypes( fetchGenotypes(srdta,1) )
m <- simpleGenotypeErrorModel(0.05, 0.01)

prob <- m["aa","Aa"]

# error model for no errors !
no_errors <- simpleGenotypeErrorModel(0,0)}

