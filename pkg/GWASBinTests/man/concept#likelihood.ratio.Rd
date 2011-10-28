\name{concept#likelihood.ratio}
\alias{concept#likelihood.ratio}
\alias{concept#likelihood.ratio}
\title{The likelihood ratio}
\description{A explained in the article, among the scores (or statistics) we consider
and use to assess the association of a bin are the likelihood ratios.}
\details{We use several different likelihood functions, that are defined on the
patients, i.e on their observed genotypes for the SNPs of the bin (cf \strong{\link{concept#bin}})  and
their different covariable values.

These likelihood functions rely upon the probability distributions of
genotypes \eqn{ P( (G)j | S, (V)m )}. Different hypotheses lead to different probability distributions
hence to different likelihood values.

To compare hypothesis, the likelihood function will be evaluated for
two different hypothesis: the H1 hypothesis considering the phenotype and
the true probability distribution \eqn{ P( (G)j | S, (V)m )}, and the null hypothesis H0
that considers that the genotypes are independent from the phenotype S:
\eqn{ P_{H_0}( (G)j | (V)m) }.

So we have two likelihoods \eqn{ L_{H_1} } and \eqn{ L_{H_0} } in function of the hypothesis/distribution considered,
that define the likelihood ratio for a given patient i:

\deqn{
LR(Patient_i) = \frac{ L_{H_1}(Patient_i) }{ L_{H_0}(Patient_i) }
}

As all patients are considered to be independently chosen, the likelihood of the set of patients available is:
\deqn{
LR( Bin_b ) = \prod_{all \ patients \ i} LR( Patient_i) )
}}
\keyword{concept}
\concept{likelihood ratio}
\alias{likelihood_ratio}

