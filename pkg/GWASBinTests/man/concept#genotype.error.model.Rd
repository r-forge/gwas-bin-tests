\name{concept#genotype.error.model}
\alias{concept#genotype.error.model}
\alias{concept#genotype.error.model}
\title{Modelization of genotyping errors}
\description{An error model is introduced with observed genotypes
\eqn{ \mathcal{O}^j } (with \eqn{ \mathcal{O}^j
\in \left\{aa,Aa,AA,\emptyset \right\} },
where \eqn{\emptyset} means that the genotype is missing):
\deqn{  P(\mathcal{O}^j | (\mathcal{G}^l)_{l \in [1,J]}) = P(\mathcal{O}^j | \mathcal{G}^j ) }}
\details{Indeed, the technology is the same for all determinations of the
same marker and there is no correlation between the genomic order
of SNPs and the localization of their probes on the genotyping
chips, so that there is no reason why the observed genotype should
depend on something else than the real genotype.

So basically, the error model for a given marker \eqn{j} is the
set of conditional probabilities

\deqn{
P( \mathcal{O}^j \backslash  \mathcal{G}^j ) =
\left(
\begin{array}{r|r|r|r}
\mathcal{O}^j \backslash  \mathcal{G}^j  & aa & Aa & AA\\
\hline

aa & & & \\
\hline
Aa & & & \\
\hline
AA & & & \\
\hline
\emptyset&  &  & \\
\end{array}
\right)
}}
\keyword{concept}
\concept{genotype error model}
\alias{genotype_error_model}

