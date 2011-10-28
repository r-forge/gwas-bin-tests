\name{concept#naive.likelihood}
\alias{concept#naive.likelihood}
\alias{concept#naive.likelihood}
\title{The Naive (aka univariate) Likelihood function}
\description{The naive or univariate likelihood is a one of the possible likelihood functions
used to assess the association of bins using likelihood ratio (cf \strong{\link{concept#likelihood.ratio}})}
\details{It is written L3 in the article, and named thereafter "naive likelihood" because it corresponds to a naive Bayesian
model. It does not model the linkage disequilibrium considering markers as independents.

For a given bin, a given set of variables (V)m, we define this likelihood function on a patient i as follows:
\deqn{
L_3(Patient_i) = \prod_{all \ markers \ j} \ \sum_{g \in {aa, Aa, AA}} p( o^j(i) | g) . p(g | (v(i))_m)
}
with \eqn{ o^j(i) } the observed genotype for marker j and patient i, and \eqn{ (v(i))_m  } the set
of variable values of patient i.

If you include the phenotype in the set of variables you end up with \eqn{L3_{H_{1}}} the likelihood under H1, and otherwise
with \eqn{L3_{H_{0}}} under H0.

As all patients are considered to be independently chosen, the likelihood of the set of patients available is:
\deqn{
L3( Bin_b ) = \prod_{all \ patients \ i} L3( Patient_i)
=  \prod_{all \ patients \ i} \ \prod_{all \ markers \ j} \ \sum_{g \in {aa, Aa, AA}} p( o^j(i) | g) . p(g | (v(i))_m)
}

We can swap the two products:
\deqn{
L3( Bin_b )
=  \prod_{all \ markers \ j} \  \prod_{all \ patients \ i} \ \sum_{g \in {aa, Aa, AA}} p( o^j(i) | g) . p(g | (v(i))_m)
}

So that if we define the likelihood L3 for a given marker j:
\deqn{
L3(Marker_j) := \prod_{all \ patients \ i} \ \sum_{g \in \{aa, Aa, AA\} } p( o^j(i) | g) . p(g | (v(i))_m
}

We can now compute likelihood L3 marker by marker, and still compute the likelihood of the bin:

\deqn{
L3( Bin_b ) =  \prod_{all \ markers \ j} \  L3( Marker_j)
}}
\keyword{concept}
\concept{naive likelihood}
\alias{univariate_likelihood,}
\alias{L3}

