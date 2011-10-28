
#library(RUnit)



test.simpleDivariateLogLikelihood <- function() {

	data("ms1")
	x_small <- ms1_gws[,998:999]
	params <- parameters(verbosity=0, nb_permutations = 0)

	cmp <- GWASBinTests:::.compareDivariateLikelihood(x_small, params)
	checkEqualsNumeric(cmp$simple_divariate_score, cmp$score_divariate)

	system.time(cmp <- GWASBinTests:::.compareDivariateLikelihood(ms1_gws, params))
	checkEqualsNumeric(cmp$simple_divariate_score, cmp$score_divariate)

	system.time(cmp <- GWASBinTests:::.compareDivariateLikelihood(ms1_gws, params, bins=ms1_bins[1:10,]))
	checkTrue( length(cmp) == 0)


	system.time(cmp <- GWASBinTests:::.compareDivariateLikelihood(ms1_gws, params, bins=ms1_bins[1:200,]))
	checkEqualsNumeric(cmp$score_divariate, cmp$simple_divariate_score)


	## using merge with different SNPS ==> missing data
	x1 <- ms1_gws[1:4, 1:2]
	x2 <- ms1_gws[5:10, 3:5]
	x <-mergeGws(x1, x2, intersected_snps_only = FALSE)

	# without covariable
	cmp <- GWASBinTests:::.compareDivariateLikelihood(x, params)
	checkEqualsNumeric(cmp$score_divariate, cmp$simple_divariate_score)

	cmp_cov <- GWASBinTests:::.compareDivariateLikelihood(x, params, covariables = data.frame(gws=x@phdata$gws) )
	checkEqualsNumeric(cmp_cov$score_divariate, cmp_cov$simple_divariate_score)

}

test.simpleDivariateLogLikelihoodForPair <- function() {
	# ===== investigate bug: it returns -Inf when e2=simpleGenotypeErrorModel(0,0)
	data("ms1")
	x_small <- ms1_gws[,998:999]
	g1 <- as.vector(fetchGenotypes(x_small, 1))
	g2 <- as.vector(fetchGenotypes(x_small, 2))
	e1 <- simpleGenotypeErrorModel(0,0.01)
	e2 <- simpleGenotypeErrorModel(0,0)
	variables <- data.frame(x_small@phdata$pop)
	llh1 <- simpleDivariateLogLikelihoodForPair(g1,e1,g2,e2, variables)
	checkTrue( ! is.infinite(llh1) )
}



test.simpleDivariatePvalue <- function() {
	## ========= TEST ==> PV~0 =====================
	g1 <-  c(rep(0L,50), rep(1L, 50))
	g2 <- rep(0L:1L,50)
	genotypes <- list(g1,g2)
	e1 <- simpleGenotypeErrorModel(0,0.01)
	e2 <- simpleGenotypeErrorModel(0,0.02)
	phenotypes <-  g1
	variables <- data.frame(phenotypes)
	em <- simpleGenotypeErrorModel(0,0)
	llh_a <- simpleDivariateLogLikelihoodForPair(g1,em,g1,em, variables)
	checkEqualsNumeric(llh_a, 0)
	llh_b <- simpleDivariateLogLikelihoodForPair(g1,e1,g2,e2, variables)
	checkTrue(llh_b < 0)

	nb_perm <- 100
	pvres <- simpleDivariatePvalue(nb_perm, genotypes, phenotypes, em)
	checkTrue(pvres$nb_better == 0)

	## ========= TEST ==> PV~1 =====================

	phenotypes <- c(rep(c(0,0,1),33),1)

	nb_perm <- 100
	pvres <- simpleDivariatePvalue(nb_perm, genotypes, phenotypes, em)
	checkTrue(pvres$nb_better > 90)

	##  =========  compare against c++ implementation  =========
	# we test H0 hypothesis : the implementations are equivalent
	# so, becasue c++ implementation is several orders of magnitude faster
	# we use it to compute a pvalue with much precision
	# then we' ll test if the simple pv follows a binomial with the c++ pvalue as
	# mean
	data("ms1")
	params <- parameters(nb_permutations=10000, verbosity=0)

	x <- ms1_gws[,40:49]

	system.time(cmp <- GWASBinTests:::.compareDivariateLikelihoodPvalues(x, params, 100))
	checkTrue( cmp$in_conf_int )

	x_small <- ms1_gws[,998:999]
	g1 <- as.vector(fetchGenotypes(x_small, 1))
	g2 <- as.vector(fetchGenotypes(x_small, 2))
	e1 <- simpleGenotypeErrorModel(0,0.01)
	e2 <- simpleGenotypeErrorModel(0,0.02)
	variables <- data.frame(x_small@phdata$pop)

	system.time(cmp <- GWASBinTests:::.compareDivariateLikelihoodPvalues(x_small, params, 1000))
	checkEqualsNumeric(cmp$score_divariate, cmp$simple_divariate_score)
	checkTrue( cmp$in_conf_int )


	# make a merge with different SNPs (and samples) and compute pvalue with study covariable
	x1 <- ms1_gws[1:100, 37:47]
	ph1 <- rep(0:1, 50)
	x2 <- ms1_gws[101:300, 101:120]
	ph2 <- rep(0:1, 100)
	x <- mergeGws(x1, x2, intersected_snps_only = FALSE)
	ph <- c(ph1, ph2)

	params <- parameters(nb_permutations=1000, verbosity=0)

	system.time(cmp_cov <- GWASBinTests:::.compareDivariateLikelihoodPvalues(x, params, 10, ph,
					covariables = data.frame(phdata(x)$gws)) )
	checkTrue( ! is.nan(cmp_cov$score_divariate) && ! is.na(cmp_cov$score_divariate)
					&& cmp_cov$pv_divariate < 1)
	checkEqualsNumeric( cmp_cov$score_divariate, cmp_cov$simple_divariate_score )
	checkTrue( cmp_cov$in_conf_int)

}

