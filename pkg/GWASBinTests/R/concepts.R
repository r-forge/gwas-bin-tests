#' A bin is a genomic region, defined by its chromosome and is start and end positions.
#' It corresponds to a contiguous set of SNPs.
#'
#' Notations:
#'
#' B is the number of bins.
#' A given bin contains J SNPs.
#'
#'
#' @name concept#bin
#' @title Definition of bins
#' @keywords concept
#' @concept bin
#' @aliases bin
NULL



#' A explained in the article, among the scores (or statistics) we consider
#' and use to assess the association of a bin are the likelihood ratios.
#'
#' We use several different likelihood functions, that are defined on the
#' patients, i.e on their observed genotypes for the SNPs of the bin (cf \strong{\link{concept#bin}})  and
#' their different covariable values.
#'
#' These likelihood functions rely upon the probability distributions of
#' genotypes \eqn{ P( (G)j | S, (V)m )}. Different hypotheses lead to different probability distributions
#' hence to different likelihood values.
#'
#' To compare hypothesis, the likelihood function will be evaluated for
#' two different hypothesis: the H1 hypothesis considering the phenotype and
#' the true probability distribution \eqn{ P( (G)j | S, (V)m )}, and the null hypothesis H0
#' that considers that the genotypes are independent from the phenotype S:
#' \eqn{ P_{H_0}( (G)j | (V)m) }.
#'
#' So we have two likelihoods \eqn{ L_{H_1} } and \eqn{ L_{H_0} } in function of the hypothesis/distribution considered,
#' that define the likelihood ratio for a given patient i:
#'
#' \deqn{
#' 	LR(Patient_i) = \frac{ L_{H_1}(Patient_i) }{ L_{H_0}(Patient_i) }
#' }
#'
#' As all patients are considered to be independently chosen, the likelihood of the set of patients available is:
#' \deqn{
#' 	LR( Bin_b ) = \prod_{all \ patients \ i} LR( Patient_i) )
#' }
#'
#'
#' @name concept#likelihood.ratio
#' @title The likelihood ratio
#' @keywords concept
#' @concept likelihood ratio
#' @aliases likelihood_ratio
NULL

#' The naive or univariate likelihood is a one of the possible likelihood functions
#' used to assess the association of bins using likelihood ratio (cf \strong{\link{concept#likelihood.ratio}})
#'
#' It is written L3 in the article, and named thereafter "naive likelihood" because it corresponds to a naive Bayesian
#' model. It does not model the linkage disequilibrium considering markers as independents.
#'
#' For a given bin, a given set of variables (V)m, we define this likelihood function on a patient i as follows:
#' \deqn{
#' L_3(Patient_i) = \prod_{all \ markers \ j} \ \sum_{g \in {aa, Aa, AA}} p( o^j(i) | g) . p(g | (v(i))_m)
#' }
#' with \eqn{ o^j(i) } the observed genotype for marker j and patient i, and \eqn{ (v(i))_m  } the set
#' of variable values of patient i.
#'
#'  If you include the phenotype in the set of variables you end up with \eqn{L3_{H_{1}}} the likelihood under H1, and otherwise
#' 	with \eqn{L3_{H_{0}}} under H0.
#'
#' As all patients are considered to be independently chosen, the likelihood of the set of patients available is:
#' \deqn{
#' 	L3( Bin_b ) = \prod_{all \ patients \ i} L3( Patient_i)
#' 	=  \prod_{all \ patients \ i} \ \prod_{all \ markers \ j} \ \sum_{g \in {aa, Aa, AA}} p( o^j(i) | g) . p(g | (v(i))_m)
#' }
#'
#' We can swap the two products:
#' \deqn{
#' 	L3( Bin_b )
#' 	=  \prod_{all \ markers \ j} \  \prod_{all \ patients \ i} \ \sum_{g \in {aa, Aa, AA}} p( o^j(i) | g) . p(g | (v(i))_m)
#' }
#'
#' So that if we define the likelihood L3 for a given marker j:
#' \deqn{
#' 	L3(Marker_j) := \prod_{all \ patients \ i} \ \sum_{g \in \{aa, Aa, AA\} } p( o^j(i) | g) . p(g | (v(i))_m
#' }
#'
#' We can now compute likelihood L3 marker by marker, and still compute the likelihood of the bin:
#'
#' \deqn{
#' 	L3( Bin_b ) =  \prod_{all \ markers \ j} \  L3( Marker_j)
#' }
#'
#' @name concept#naive.likelihood
#' @title The Naive (aka univariate) Likelihood function
#' @keywords concept
#' @concept naive likelihood
#' @aliases univariate_likelihood, L3
NULL

#' The  two-marker sliding windows (aka divariate) likelihood is a one of the possible likelihood functions
#' used to assess the association of bins using likelihood ratio (cf \strong{\link{concept#likelihood.ratio}})
#'
#' It is written L2 in the article.
#'
#' For a given bin, a given set of variables (V)m, we define this likelihood function on a patient i as follows:
#' \deqn{
#' L_2(Patient_i) = \prod_{  2 \le j \le J } \ \sum_{(g_1, g_2) \in {[aa, Aa, AA]}^2}
#'    p( o^j(i) | g_2) . p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m) . p( o^{j-1}(i) | g_1)
#' }
#' with \eqn{ o^j(i) } the observed genotype for marker j and patient i, \eqn{ (v(i))_m  } the set
#' of variable values of patient i and \eqn{p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m)} the probability of
#' having the two genotypes for markers j-1 and j.
#'
#' If you include the phenotype in the set of variables you end up with \eqn{L2_{H_{1}}} the likelihood under H1, and otherwise
#' 	with \eqn{L2_{H_{0}}} under H0.
#'
#' As all patients are considered to be independently chosen, the likelihood of the set of patients available is:
#' \deqn{
#' 	L2( Bin_b ) = \prod_{all \ patients \ i} L2( Patient_i)
#' 	=  \prod_{all \ patients \ i} \ \prod_{2 \le j \le J }
#' \ \sum_{(g_1, g_2) \in {[aa, Aa, AA]}^2}
#'    p( o^j(i) | g_2) . p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m) . p( o^{j-1}(i) | g_1)
#' }
#'
#' We can swap the two products:
#' \deqn{
#' 	L2( Bin_b )
#' 	=  \prod_{2 \le j \le J } \  \prod_{all \ patients \ i} \ \sum_{(g_1, g_2) \in {[aa, Aa, AA]}^2}
#'    p( o^j(i) | g_2) . p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m) . p( o^{j-1}(i) | g_1)
#' }
#'
#'
#' So that if we define the likelihood L2 for a pair of consecutive markers (j-1, j):
#' \deqn{
#' 	L2(j-1,j) := \prod_{all \ patients \ i} \ \sum_{(g_1, g_2) \in {[aa, Aa, AA]}^2}   p( o^j(i) | g_2) . p(G_{(j-1, j)}=(g_1, g_2) | (v(i))_m) . p( o^{j-1}(i) | g_1)
#' }
#'
#' We can now compute likelihood L2 successively for each pair of consecutive markers, and still compute the likelihood of the bin:
#'
#' \deqn{
#' 	L2( Bin_b ) =  \prod_{2 \le j \le J } \  L2(j-1,j)
#' }
#'
#' @name concept#divariate.likelihood
#' @title The two-marker sliding windows (aka divariate) Likelihood function
#' @keywords concept
#' @concept divariate likelihood
#' @aliases divariate_likelihood, L2
NULL

#' The allelic Sum Of Scores statistic.
#'
#' It is based on the Pearson test-statistic on a contingency table (a.k.a chi-squared test) \eqn{X^2}.
#'
#' For a bin b, a marker j, a given set of variables values \eqn{(v)_m} selects a subset of the individuals I.
#' On this subset I, we note the contingency table cell values counting the number of alleles per phenotype value
#' \eqn{O(a,p)} where \eqn{O(a,p)} is the number of allele a for phenotype p.
#'
#' Then the allelic score (Pearson score) for this table is:
#' \deqn{
#' X_{all}^2(b,j,(v)_m) = \sum_{ \{ a \in \{a,A\},p \in P,  E(a,p) \ne 0 \} } \ \frac{(O(a,p)-E(a,p))^2}{E(a,p)}
#' }
#'
#' so that
#'
#' \deqn{
#' X_{all}^2(Bin\ b)=\sum_{marker_j \in b} \sum_{(v)_m \in (V)_m} X_{all}^2(b,j,(v)_m)
#' }
#'
#'
#'
#' @name concept#allelic.sum.score
#' @title The allelic Sum of scores bin statistic
#' @keywords concept
#' @concept allelic SumScore
#' @aliases allelic_sum_score
NULL

#' The genotypic Sum Of Scores statistic.
#'
#' It is the equivalent of the \link{concept#allelic.sum.score} but computed on genotypic tables
#' instead of allelic ones.
#'
#' It is based on the Pearson test-statistic on a contingency table (a.k.a chi-squared test) \eqn{X^2}.
#'
#' For a bin b, a marker j, a given set of variables values \eqn{(v)_m} selects a subset of the individuals I.
#' On this subset I, we note the contingency table cell values counting the number of genotypes per phenotype value
#' \eqn{O(g,p)} where \eqn{O(g,p)} is the number of genotypes g for phenotype p.
#'
#' Then the genotypic score (Pearson score) for this table is:
#' \deqn{
#' X_{gen}^2(b,j,(v)_m) = \sum_{ \{ g \in \{aa,Aa, AA\},p \in P,  E(g,p) \ne 0 \} } \ \frac{(O(g,p)-E(g,p))^2}{E(g,p)}
#' }
#'
#' so that
#'
#' \deqn{
#' X_{gen}^2(Bin\ b)=\sum_{marker_j \in b} \sum_{(v)_m \in (V)_m} X_{gen}^2(b,j,(v)_m)
#' }
#'
#'
#'
#' @name concept#genotypic.sum.score
#' @title The genotypic Sum of scores bin statistic
#' @keywords concept
#' @concept genotypic SumScore
#' @aliases genotypic_sum_score
NULL

#' An error model is introduced with observed genotypes
#' \eqn{ \mathcal{O}^j } (with \eqn{ \mathcal{O}^j
#' 						\in \left\{aa,Aa,AA,\emptyset \right\} },
#' where \eqn{\emptyset} means that the genotype is missing):
#' \deqn{  P(\mathcal{O}^j | (\mathcal{G}^l)_{l \in [1,J]}) = P(\mathcal{O}^j | \mathcal{G}^j ) }
#'
#' Indeed, the technology is the same for all determinations of the
#' same marker and there is no correlation between the genomic order
#' of SNPs and the localization of their probes on the genotyping
#' chips, so that there is no reason why the observed genotype should
#' depend on something else than the real genotype.
#'
#' So basically, the error model for a given marker \eqn{j} is the
#' set of conditional probabilities
#'
#' \deqn{
#'  P( \mathcal{O}^j \backslash  \mathcal{G}^j ) =
#'  \left(
#'  \begin{array}{r|r|r|r}
#'  \mathcal{O}^j \backslash  \mathcal{G}^j  & aa & Aa & AA\\
#'  \hline
#'
#'  aa & & & \\
#'  \hline
#'  Aa & & & \\
#'  \hline
#'  AA & & & \\
#'  \hline
#'  \emptyset&  &  & \\
#'  \end{array}
#'  \right)
#' }
#'
#'
#' @name concept#genotype.error.model
#' @title Modelization of genotyping errors
#' @keywords concept
#' @concept genotype error model
#' @aliases genotype_error_model

NULL

#' To obtain more regular estimates, a constant is added to all cell counts.
#' It is a Dirichlet prior on parameters. This constant is
#' can be chosen to be \eqn{C=\alpha . \bar{n}}, where
#' \eqn{\alpha} is the chosen error rate and
#' \eqn{\bar{n}} is the mean number of individuals
#' per cell. This constant means that uncertainty on low cell counts
#' is high, not only because of randomness, but also because of
#' genotyping errors.
#'
#' @name concept#regularization.of.table.estimates
#' @title Regularizations of the estimates from contingency tables
#' @keywords concept
#' @concept regularization of table estimates
#' @aliases regularization_of_table_estimates

NULL