\name{concept#divariate.likelihood}
\alias{concept#divariate.likelihood}
\alias{concept#divariate.likelihood}
\title{The two-marker sliding windows (aka divariate) Likelihood function}
\description{The  two-marker sliding windows (aka divariate) likelihood is a one of the possible likelihood functions
used to assess the association of bins using likelihood ratio (cf \strong{\link{concept#likelihood.ratio}})}
\details{It is written L2 in the article.

For a given bin, a given set of variables (V)m, we define this likelihood function on a patient i as follows:
\deqn{
L_2(Patient_i) = \prod_{  2 \le j \le J } \ \sum_{(g_1, g_2) \in {[aa, Aa, AA]}^2}
p( o^j(i) | g_2) . p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m) . p( o^{j-1}(i) | g_1)
}
with \eqn{ o^j(i) } the observed genotype for marker j and patient i, \eqn{ (v(i))_m  } the set
of variable values of patient i and \eqn{p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m)} the probability of
having the two genotypes for markers j-1 and j.

If you include the phenotype in the set of variables you end up with \eqn{L2_{H_{1}}} the likelihood under H1, and otherwise
with \eqn{L2_{H_{0}}} under H0.

As all patients are considered to be independently chosen, the likelihood of the set of patients available is:
\deqn{
L2( Bin_b ) = \prod_{all \ patients \ i} L2( Patient_i)
=  \prod_{all \ patients \ i} \ \prod_{2 \le j \le J }
\ \sum_{(g_1, g_2) \in {[aa, Aa, AA]}^2}
p( o^j(i) | g_2) . p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m) . p( o^{j-1}(i) | g_1)
}

We can swap the two products:
\deqn{
L2( Bin_b )
=  \prod_{2 \le j \le J } \  \prod_{all \ patients \ i} \ \sum_{(g_1, g_2) \in {[aa, Aa, AA]}^2}
p( o^j(i) | g_2) . p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m) . p( o^{j-1}(i) | g_1)
}


So that if we define the likelihood L2 for a pair of consecutive markers (j-1, j):
\deqn{
L2(j-1,j) := \prod_{all \ patients \ i} \ \sum_{(g_1, g_2) \in {[aa, Aa, AA]}^2}   p( o^j(i) | g_2) . p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m) . p( o^{j-1}(i) | g_1)
}

We can now compute likelihood L2 successively for each pair of consecutive markers, and still compute the likelihood of the bin:

\deqn{
L2( Bin_b ) =  \prod_{2 \le j \le J } \  L2(j-1,j)
}}
\keyword{concept}
\concept{divariate likelihood}
\alias{divariate_likelihood,}
\alias{L2}

