\name{concept#genotypic.sum.score}
\alias{concept#genotypic.sum.score}
\alias{concept#genotypic.sum.score}
\title{The genotypic Sum of scores bin statistic}
\description{The genotypic Sum Of Scores statistic.}
\details{It is the equivalent of the \link{concept#allelic.sum.score} but computed on genotypic tables
instead of allelic ones.

It is based on the Pearson test-statistic on a contingency table (a.k.a chi-squared test) \eqn{X^2}.

For a bin b, a marker j, a given set of variables values \eqn{(v)_m} selects a subset of the individuals I.
On this subset I, we note the contingency table cell values counting the number of genotypes per phenotype value
\eqn{O(g,p)} where \eqn{O(g,p)} is the number of genotypes g for phenotype p.

Then the genotypic score (Pearson score) for this table is:
\deqn{
X_{gen}^2(b,j,(v)_m) = \sum_{ \{ g \in \{aa,Aa, AA\},p \in P,  E(g,p) \ne 0 \} } \ \frac{(O(g,p)-E(g,p))^2}{E(g,p)}
}

so that

\deqn{
X_{gen}^2(Bin\ b)=\sum_{marker_j \in b} \sum_{(v)_m \in (V)_m} X_{gen}^2(b,j,(v)_m)
}}
\keyword{concept}
\concept{genotypic SumScore}
\alias{genotypic_sum_score}

