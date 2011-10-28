
#library(RUnit)



test.simpleAllelicSumScores <- function() {
	g <-  rep(0L,100)
	phenotypes <-  rep(0:1, 50)
	s <- simpleAllelicSumScores(list(g), phenotypes)
	checkEqualsNumeric(0, s)

	g <-  rep(0L:1L, 50)
	phenotypes <- g
	s <- simpleAllelicSumScores(list(g), phenotypes)
	# checked with :  chisq.test
	checkEqualsNumeric(s, chisq.test(matrix(c(100, 0, 50, 50), byrow = T, nrow=2), correct=F)$statistic)

	# all missing ==> NA
	g <-  rep(-1L,100)
	phenotypes <-  rep(0:1, 50)
	s <- simpleAllelicSumScores(list(g), phenotypes)
	checkTrue( s == 0)

	## exactly average
	g <-  rep(-1L:2L, 25)
	phenotypes <- rep(0L:1L,50)
	s <- simpleAllelicSumScores(list(g), phenotypes)
	checkEqualsNumeric(0, s)

	## ==== check vs nnbc++ ====
	data("ms1")
	params <- parameters(verbosity=0, nb_permutations = 0)

	# one SNP
	x1 <- ms1_gws[,333]
	cmp <- GWASBinTests:::.compareAllelicSumScores(x1, params)
	checkEqualsNumeric(cmp$score_allelic, cmp$simple_allelic_score)

	# 2 SNPs
	x_small <- ms1_gws[,998:999]
	cmp <- GWASBinTests:::.compareAllelicSumScores(x_small, params)
	checkEqualsNumeric(cmp$score_allelic, cmp$simple_allelic_score)

	# 100 SNPs
	x100 <- ms1_gws[,557:656]
	cmp <- GWASBinTests:::.compareAllelicSumScores(x100, params)
	checkEqualsNumeric(cmp$score_allelic, cmp$simple_allelic_score)

	# on 51 bins
	system.time(cmp <- GWASBinTests:::.compareAllelicSumScores(ms1_gws, params, bins=ms1_bins[1:200,]))
	checkEqualsNumeric(cmp$score_allelic, cmp$simple_allelic_score)

	# take one SNP and only retain the NA genotypes
	only_missing <- ms1_gws[,2]
	only_missing <- only_missing[is.na(as.numeric(only_missing)),]
	system.time(cmp <- GWASBinTests:::.compareAllelicSumScores(only_missing, params) )


}



test.simpleAllelicSumScoresPvalue <- function() {
	## ==== check vs nnbc++ ====
	data("ms1")
	params <- parameters(verbosity=0, nb_permutations = 10000)

	# one SNP
	x1 <- ms1_gws[,333]
	system.time( cmp <- GWASBinTests:::.compareAllelicSumScoresPvalues(x1, params, 1000) )
	checkTrue(cmp$in_conf_int)

	# on bin 68
	bin68 <- ms1_bins[68,]
	system.time( cmp <- GWASBinTests:::.compareAllelicSumScoresPvalues(ms1_gws, params, 200, bins=bin68) )
	checkTrue(cmp$in_conf_int)


	# on some bins
	system.time(cmp <- GWASBinTests:::.compareAllelicSumScoresPvalues(ms1_gws, params, 100, bins=ms1_bins[1:75,]))
	checkEqualsNumeric(cmp$score_allelic, cmp$simple_allelic_score)
	checkTrue( all(cmp$in_conf_int) )

	# make a merge with different SNPs (and samples) and compute pvalue with study covariable
	x1 <- ms1_gws[1:100, 37:47]
	ph1 <- rep(0:1, 50)
	x2 <- ms1_gws[101:300, 101:120]
	ph2 <- rep(0:1, 100)
	x <- mergeGws(x1, x2, intersected_snps_only = FALSE)
	ph <- c(ph1, ph2)

	params <- parameters(nb_permutations=1000, verbosity=0)

	system.time(cmp_cov <- GWASBinTests:::.compareAllelicSumScoresPvalues(x, params, 10, ph,
					covariables = data.frame(phdata(x)$gws)) )
	checkTrue( ! is.nan(cmp_cov$score_allelic) && ! is.na(cmp_cov$simple_allelic_score)
					&& cmp_cov$pv_allelic < 1 )
	checkEqualsNumeric( cmp_cov$score_allelic, cmp_cov$simple_allelic_score )
	checkTrue( cmp_cov$in_conf_int)

}

