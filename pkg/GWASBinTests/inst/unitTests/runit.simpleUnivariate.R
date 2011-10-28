
#library(RUnit)



test.simpleUnivariateLogLikelihood <- function() {

	## ========= TEST 1 =====================

	genotypes1 <- as.integer(c(0, 1))
	genotypes2 <- as.integer( c(2, 2 ))
	genotypes_list <- list(genotypes1, genotypes2)

	e1 <- simpleGenotypeErrorModel(0.05, 0.01)
	e2 <- simpleGenotypeErrorModel(0,0)
	error_models <- list(e1, e2)

	phenotypes <-  data.frame(phenos=c( 0, 1 ))
	llh <- simpleUnivariateLogLikelihood(genotypes_list, phenotypes, error_models,  regCoeff=0.05)

	# check it manually
	llh1 <- simpleUnivariateLogLikelihoodForMarker(genotypes1, phenotypes, e1, regCoeff=0.05)
	llh2 <- simpleUnivariateLogLikelihoodForMarker(genotypes2, phenotypes, e2, regCoeff=0.05)
	checkEqualsNumeric(llh, llh1 + llh2)

	# calling without variables
	llh0 <- simpleUnivariateLogLikelihood(genotypes_list, data.frame(), error_models,  regCoeff=0.05)

	## ========= TEST 2 =====================

	genotypes1 <- as.integer(c(0, 1))
	genotypes2 <- as.integer( c(2, 2 ))
	phenotypes <-  c( 0, 1 )
	h1 <- data.frame(h1=phenotypes)
	genos <- list(genotypes1, genotypes2)
	error_model <- matrix(c(0.95, 0.05, 0
					, 0.05, 0.90, 0.05
					,0, 0.05, 0.95
					,0, 0, 0 ),nrow=4, ncol=3, byrow = TRUE)

	expected <- log(prod(c(1*0.95*1*0.9, 1*0.95*1*0.95)))
	llh <- simpleUnivariateLogLikelihood(genos, h1, list(error_model, error_model))
	checkEqualsNumeric(llh, expected)

	# test alternate syntax for unique error model
	llhbis <- simpleUnivariateLogLikelihood(genos, h1, error_model)
	checkEqualsNumeric(llh, llhbis)

	# ==== TEST against c++ implementation =====
	data("ms1")

	params <- parameters(nb_permutations=1, verbosity=0)

	# 100 SNPs
	x <- ms1_gws[,23:122]
	cmp <- GWASBinTests:::.compareUnivariateLikelihood(x, params)
	checkEqualsNumeric(cmp$score_univariate, cmp$simple_univariate_score)

	# ==== TEST merging of identical SNPs ====
	x1 <- ms1_gws[,10]
	x2 <- x1
	# rename x2 samples
	x2ids <- paste("dummy", idnames(x2), sep="_")
	x2@phdata$id <- x2ids
	x2@gtdata@idnames <- x2ids

	x <-mergeGws(x1, x2)
	cmp <- GWASBinTests:::.compareUnivariateLikelihood(x, params)
	cmp1 <- GWASBinTests:::.compareUnivariateLikelihood(x1, params)

	checkEqualsNumeric(cmp$score_univariate, cmp$simple_univariate_score)
	checkEqualsNumeric(cmp1$score_univariate, cmp1$simple_univariate_score)
	checkEqualsNumeric(cmp$simple_univariate_score, 2*cmp1$simple_univariate_score)


# ==== TEST merging of different SNPs ====

	x1 <- ms1_gws[,10]
	x2 <- ms1_gws[,189]
	# rename x2 snp
	x2bis <- x2
	x2ids <- paste("dummy", idnames(x2), sep="_")
	x2bis@phdata$id <- x2ids
	x2bis@gtdata@idnames <- x2ids
	x2bis@gtdata@snpnames <- snpnames(x1)
	x2bis@gtdata@chromosome <- factor(chromosome(x1))
	x2bis@gtdata@coding <- x1@gtdata@coding
	x2bis@gtdata@map <- x1@gtdata@map

	x <-mergeGws(x1, x2bis)

	cmp <- GWASBinTests:::.compareUnivariateLikelihood(x, params)
	cmp_cov <- GWASBinTests:::.compareUnivariateLikelihood(x, params, covariables = data.frame(gws=x@phdata$gws) )
	cmp1 <- GWASBinTests:::.compareUnivariateLikelihood(x1, params)
	cmp2 <- GWASBinTests:::.compareUnivariateLikelihood(x2bis, params)

	checkEqualsNumeric(cmp_cov$score_univariate, cmp_cov$simple_univariate_score)

	checkEqualsNumeric(cmp$score_univariate, cmp$simple_univariate_score)
	checkEqualsNumeric(cmp1$score_univariate, cmp1$simple_univariate_score)
	checkEqualsNumeric(cmp2$score_univariate, cmp2$simple_univariate_score)


	## !!! we can see here that the likelihood of the merged is not simply
	# the sum of the likelihoods !
	checkException( checkEqualsNumeric(cmp$simple_univariate_score, cmp1$simple_univariate_score+cmp2$simple_univariate_score), silent = TRUE )


	## using merge with different SNPS ==> missing data
	x1 <- ms1_gws[1:4, 1]
	x2 <- ms1_gws[5:10, 2]
	x <-mergeGws(x1, x2, intersected_snps_only = FALSE)

	# without covariable
	cmp <- GWASBinTests:::.compareUnivariateLikelihood(x, params)
	checkEqualsNumeric(cmp$score_univariate, cmp$simple_univariate_score)

	cmp_cov <- GWASBinTests:::.compareUnivariateLikelihood(x, params, covariables = data.frame(gws=x@phdata$gws) )
	checkEqualsNumeric(cmp_cov$score_univariate, cmp_cov$simple_univariate_score)
}

test.simpleUnivariateLogLikelihoodForMarker <- function() {
	error_model <- matrix(c(0.95, 0.05, 0
							, 0.05, 0.90, 0.05
							,0, 0.05, 0.95
							,0, 0, 0 ),nrow=4, ncol=3, byrow = TRUE)
	genotypes1 <- as.integer(c(0, 1))
	genotypes2 <- as.integer( c(2, 2 ) )
	phenotypes <-  c( 0, 1 )
	h1 <- data.frame(h1=phenotypes)
	h0 <- data.frame(h0=rep(0,2))

	hyp <- list(h1,h0)
	genos <- list(genotypes1, genotypes2)

	expected <- c(1*0.95*1*0.9, (0.95+0.05)/2 * (0.9+0.05)/2,
			1*0.95*1*0.95, 1*0.95 * 1*0.95)

	llhs <- lapply(genos, function(g) lapply(hyp, function(h) {
		llh <- simpleUnivariateLogLikelihoodForMarker(g, h, error_model)
	}  ))
	lhs <- exp(unlist(llhs))
	checkEqualsNumeric(lhs, expected)

	## test when all genotypes are missing
	## in that case the likelihood should not be defined: NaN
	## test when all genotypes are missing
	## in that case the likelihood should not be defined: NaN

	genotypes <- as.integer( c(-1, -1) )
	h1 <- data.frame(h1=phenotypes)
	llh <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, error_model)
	checkTrue( llh == 0 )


	## test caricatural example with null error model
	## and genotypes=(0x50,1x50)

	null_model <- simpleGenotypeErrorModel(0,0)
	genotypes <- c(rep(0L,50), rep(1L,50) )
	h0 <- data.frame(h0=rep(0,100))
	h1 <- data.frame(h1=genotypes)
	llh0 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h0, null_model)
	# llh0=100*ln( 1/2 )
	checkEqualsNumeric(llh0, 100*log(0.5))
	llh1 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, null_model)
	# llh1=100*ln(1)=0
	checkEqualsNumeric(llh1, 0)

	h2 <- rev(h1)
	llh2 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h2, null_model)
	checkEqualsNumeric(llh2, llh1)

	h3 <- data.frame(pheno=rep(c(0,1), 50))
	llh3 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h3, null_model)
	# llh3 = 4 * 25*ln( 1/2 ) == llh0
	checkEqualsNumeric(llh3, llh0)

	# no variables = one constant variable == h0
	llh_novar <- simpleUnivariateLogLikelihoodForMarker(genotypes, data.frame(), null_model)
	checkEqualsNumeric(llh_novar, llh0)


	# all missing data
	genotypes <- c(rep(3L, 1000) )
	h1 <- data.frame(pheno=rep(c(0,1), 500))
	llh0 <- simpleUnivariateLogLikelihoodForMarker(genotypes, data.frame(), null_model)
	checkEqualsNumeric(llh0, 0)
	llh1 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, null_model)
	checkEqualsNumeric(llh1, 0)

	# half missing data
	genotypes <- rep(0L:3L, 25)
	h1 <- data.frame(pheno=rep(c(0,1), 50))
	llh0 <- simpleUnivariateLogLikelihoodForMarker(genotypes, data.frame(), null_model)
	checkEqualsNumeric(llh0, 75*log( 1/3 ))
	llh1 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, null_model)
	checkEqualsNumeric(llh1, 50*log(1/2))

	# real example
	genotypes = as.integer(c(rep(0,20), rep(1,20), rep(2, 1), rep(3,3), rep(0,43), rep(1,10), rep(2,3)  ))
	h1 <- data.frame(pheno=c( rep(0, 44), rep(1, 56)) )
	llh0 <- simpleUnivariateLogLikelihoodForMarker(genotypes, data.frame(), null_model)
	llh1 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, null_model)
	checkEqualsNumeric(llh0, 63*log(63/97)+30*log(30/97)+4*log(4/97))
	checkEqualsNumeric(llh1, 2*20*log(20/41) + log(1/41) + 43*log(43/56) + 10*log(10/56) + 3*log(3/56) )
	checkTrue( llh1 > llh0 )

	# influence of regCoeff
	llh0_05 <- simpleUnivariateLogLikelihoodForMarker(genotypes, data.frame(), null_model, regCoeff = 0.05)
	tot <- 97*1.05
	cc <- 0.05*97/3
	checkEqualsNumeric(llh0_05, 63*log((63+cc)/tot)+30*log((30+cc)/tot)+4*log((4+cc)/tot))

	llh1_05 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, null_model, regCoeff = 0.05)
	checkTrue( llh1_05 != llh1)
	llh1_1 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, null_model, regCoeff = 0.1)
	checkTrue( llh1_1 != llh1)
	llh1_5 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, null_model, regCoeff = 0.5)
	checkTrue( llh1_5 != llh1)

	# influence of the error model
	llh0_bis <- simpleUnivariateLogLikelihoodForMarker(genotypes, data.frame(), simpleGenotypeErrorModel(0.05, 0.05))
	checkTrue( llh0_bis != llh0)

	## test with co-variables
	# try something extreme, all genotype 0 are women, others are men
	gender <- (genotypes == 0) + 1
	h0 <- data.frame(gender=gender)
	h1 <- data.frame(pheno=c( rep(0, 44), rep(1, 56)), gender=gender )

	llh0 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h0, null_model)
	checkEqualsNumeric(llh0, 30*log(30/34) + 4*log(4/34 ))
	llh1 <- simpleUnivariateLogLikelihoodForMarker(genotypes, h1, null_model)
	checkEqualsNumeric(llh1, 20*log(20/21) + log(1/21 ) + 10*log(10/13) + 3*log(3/13 ))
}



test.simpleUnivariateLogLikelihoodRatio <- function() {

	error_model <- matrix(c(0.95, 0.05, 0
					, 0.05, 0.90, 0.05
					,0, 0.05, 0.95
					,0, 0, 0 ),nrow=4, ncol=3, byrow = TRUE)

	g1 <- as.integer(c(0, 1))
	g2 <- as.integer(c( 2, 2 ))

	genotypes_list <- list(g1,g2)
	phenotypes <-  c( 0, 1 )
	llr <- simpleUnivariateLogLikelihoodRatio(genotypes_list, phenotypes, error_model)
	checkEqualsNumeric(exp(llr), 3.6)

	# check with different regCoeff
	llr2 <- simpleUnivariateLogLikelihoodRatio(genotypes_list, phenotypes, error_model, regCoeff = 0.05)
	checkTrue( ! isTRUE(all.equal(llr, llr2)) )

	# =========  check basic multi-study with one snp duplicated =============
	# merge on same snp with distinct ids
	# theoretically, the llrx1==llrx2 and llrx== 2*llrx1
	data(srdta) # genabel dataset

	srdta@phdata$pop <- srdta@phdata$bt # we use the pop variable as phenotype
	gws <- asGws(srdta, assignPopNATo = 0)
	x1 <- gws[,2]
	x2 <- gws[,2]
	# change ids
	x1@gtdata@idnames <- paste("x1", x1@gtdata@idnames, sep="")
	x1@phdata$id <- x1@gtdata@idnames
	x <- mergeGws(x1, x2)

	# so x1 and x2 are basically the same dataset, and x is the concatenation of the 2

	m <- simpleGenotypeErrorModel(0.05,0.01)
	g <- fetchGenotypesAsList(x)
	g1 <- fetchGenotypesAsList(x1)
	g2 <- fetchGenotypesAsList(x2)
	p <- as.vector(unlist(phdata(x)["pop"]))
	p1 <- as.vector(unlist(phdata(x1)["pop"]))
	p2 <- as.vector(unlist(phdata(x1)["pop"]))

	llrx <- simpleUnivariateLogLikelihoodRatio(g,p , m)
	llrx1 <- simpleUnivariateLogLikelihoodRatio(g1, p1, m)
	llrx2 <- simpleUnivariateLogLikelihoodRatio(g2, p2, m)
	checkEquals(llrx1, llrx2)
	checkEquals(llrx, 2*llrx1)

	# use the study/collection as covariable ==> should give the same
	llrx_by_study <- simpleUnivariateLogLikelihoodRatio(g, p, m, covariables = data.frame(study=x@phdata$gws))
	checkEquals(llrx, llrx_by_study)


	## check co-variable impact: compare with and w/o sex
	srdta@phdata$pop <- srdta@phdata$bt
	gws <- asGws(srdta, assignPopNATo = 0)
	x <- gws[,2]
	m <- simpleGenotypeErrorModel(0.05,0.01)
	phenos <- phdata(x)["pop"]
	llrx_nosex <- simpleUnivariateLogLikelihoodRatio(fetchGenotypesAsList(x), phenos, m)
	llrx_sex <- simpleUnivariateLogLikelihoodRatio(fetchGenotypesAsList(x), phenos, m, covariables = phdata(x)["sex"] )
	checkTrue( llrx_sex >= llrx_nosex)

	## try random co-variable: any co-variable will always
	## increase the likelihood ratio
	g <- fetchGenotypesAsList(x)
	llrx_norandom <- simpleUnivariateLogLikelihoodRatio(g, phenos, m)
#
#	llrx_randoms <- replicate(20,
#			simpleUnivariateLogLikelihoodRatio(g, phenos, m, covariables = data.frame(random=as.integer(runif(nids(x), 0,1)+0.5)))
#	)
#	checkTrue( all( llrx_randoms >= llrx_norandom) )

	## try constant covariable

	llrx_constant <- simpleUnivariateLogLikelihoodRatio(g, phenos, m, covariables = data.frame(constant=rep(0,2500)) )
	checkEqualsNumeric(llrx_constant, llrx_norandom)

	# =========  compare against  univariateLikelihoodRatio=============
#	#bin <- list(chr_name=1, start=816845, end=859024)
#	m <- simpleGenotypeErrorModel(0.05,0.01)
#	genotypes <- fetchGenotypes(srdta, 1:10)
#	genotypes_list <- fetchGenotypesAsList(srdta, 1:10)
#	vars <- phdata(srdta)[c("bt", "sex")]
#	sex <- phdata(srdta)["sex"]
#	phenos <-  phdata(srdta)["bt"]
#	llr1 <- simpleUnivariateLogLikelihoodRatio(genotypes_list, phenos,  m, covariables = sex)
#	llr2 <- univariateLikelihoodRatio(genotypes ,  vars, m)
#	checkEquals(llr1, llr2)


	## ======= merge on same samples, different SNPs =======
	## not allowed for now via mergeGws
	## Should be the sum, works fine for univariate
	x1 <- gws[1:100,1:10]
	x2 <- gws[1:100, 20:25]
	m <- simpleGenotypeErrorModel(0.05,0.01)
	x <- merge(x1, x2)
	x@phdata$pop = x1@phdata$pop
	llrx <- simpleUnivariateLogLikelihoodRatio(fetchGenotypesAsList(x), phdata(x)["pop"], m)
	llrx1 <- simpleUnivariateLogLikelihoodRatio(fetchGenotypesAsList(x1), phdata(x1)["pop"], m)
	llrx2 <- simpleUnivariateLogLikelihoodRatio(fetchGenotypesAsList(x2), phdata(x2)["pop"], m)
	checkEquals(llrx, llrx1+llrx2)

}


test.simpleUnivariatePvalue <- function() {
	## ========= TEST ==> PV~0 =====================
	genotypes <- list( c(rep(0L,50), rep(1L, 50)) )
	em <- simpleGenotypeErrorModel(0,0)
	phenotypes <-  genotypes[[1]]
	llh <- simpleUnivariateLogLikelihood(genotypes, data.frame(phenotypes), em)
	checkEqualsNumeric(llh, 0)

	nb_perm <- 100
	pvres <- simpleUnivariatePvalue(nb_perm, genotypes, phenotypes, em)
	checkTrue(pvres$nb_better <= 1)

	## ========= TEST ==> PV~1 =====================
	genotypes <- list( rep(c(0L,1L), 50) )
	em <- simpleGenotypeErrorModel(0,0)
	phenotypes <-  c(rep(0L,50), rep(1L, 50))
	llh <- simpleUnivariateLogLikelihood(genotypes, data.frame(phenotypes), em)
	pvres <- simpleUnivariatePvalue(100, genotypes, phenotypes, em)
	checkTrue(pvres$nb_better >= 99)


	## compare against c++ implementation
	# we test H0 hypothesis : the implementations are equivalent
	# so, becasue c++ implementation is several orders of magnitude faster
	# we use it to compute a pvalue with much precision
	# then we' ll test if the simple pv follows a binomial with the c++ pvalue as
	# mean
	data("ms1")
	params <- parameters(nb_permutations=10000, verbosity=0)

	x <- ms1_gws[,40:49]
	system.time(cmp <- GWASBinTests:::.compareUnivariateLikelihood(x, params))
	checkEqualsNumeric(cmp$score_univariate, cmp$simple_univariate_score)

	system.time(cmp <- GWASBinTests:::.compareUnivariateLikelihoodPvalues(x, params, 100))
	checkTrue( cmp$in_conf_int )

	# make a merge with different SNPs (and samples) and compute pvalue with study covariable
	x1 <- ms1_gws[1:100, 37:47]
	ph1 <- rep(0:1, 50)
	x2 <- ms1_gws[101:300, 101:120]
	ph2 <- rep(0:1, 100)
	x <- mergeGws(x1, x2, intersected_snps_only = FALSE)
	ph <- c(ph1, ph2)

	params <- parameters(nb_permutations=1000, verbosity=0)

	system.time(cmp_cov <- GWASBinTests:::.compareUnivariateLikelihoodPvalues(x, params, 100, ph,
					covariables = data.frame(phdata(x)$gws)) )
	checkTrue( ! is.nan(cmp_cov$score_univariate) && ! is.na(cmp_cov$score_univariate)
					&& cmp_cov$pv_univariate < 1)
	checkEqualsNumeric( cmp_cov$score_univariate, cmp_cov$score_univariate )
	checkTrue( cmp_cov$in_conf_int)

}

