\name{concept#allelic.sum.score}
\alias{concept#allelic.sum.score}
\alias{concept#allelic.sum.score}
\title{The allelic Sum of scores bin statistic}
\description{The allelic Sum Of Scores statistic.}
\details{It is based on the Pearson test-statistic on a contingency table (a.k.a chi-squared test) \eqn{X^2}.

For a bin b, a marker j, a given set of variables values \eqn{(v)_m} selects a subset of the individuals I.
On this subset I, we note the contingency table cell values counting the number of alleles per phenotype value
\eqn{O(a,p)} where \eqn{O(a,p)} is the number of allele a for phenotype p.

Then the allelic score (Pearson score) for this table is:
\deqn{
X_{all}^2(b,j,(v)_m) = \sum_{ \{ a \in \{a,A\},p \in P,  E(a,p) \ne 0 \} } \ \frac{(O(a,p)-E(a,p))^2}{E(a,p)}
}

so that

\deqn{
X_{all}^2(Bin\ b)=\sum_{marker_j \in b} \sum_{(v)_m \in (V)_m} X_{all}^2(b,j,(v)_m)
}}
\keyword{concept}
\concept{allelic SumScore}
\alias{allelic_sum_score}

