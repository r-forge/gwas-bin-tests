#library(RUnit)
# library(GWASBinTests, lib.loc = "..truetest/")
# unloadNamespace("GWASBinTests")
# library(R.utils)
# reassignInPackage("processGws", "GWASBinTests", processGws)
# rm(processGws)


test.processGws.merged.dataset <- function() {
	data(srdta)
	binsfile <-  "extdata/ms1.bins"
	bins <- readBins(binsfile)

	# restrict to the first 5 bins
	bins <- bins[1:5,]

	gws <- asGws(srdta, assignPopNATo = 0)
	# split in two
	half <- as.integer(gws@gtdata@nids/2)
	gws1 <- gws[1:half]
	gws2 <- gws[(half+1):gws@gtdata@nids]

	# set the pop phenotype to 0 for first and 1 for second
	gws1@phdata$pop <- 1
	gws2@phdata$pop <- 2


	# merge it
	merged <- mergeGws(gws1, gws2)

	# run nnbc on both
	params <- parameters(
			seed = 0,
			verbosity = 0,
			nb_permutations = 100)

	res1 <- processGws(gws1, bins, params )
	res2 <- processGws(gws2, bins, params )
	res.merged <- processGws(merged, bins, params)

	# gws1 and 2 have only one phenotype ==> all pvalues == 1
	checkEquals(res1[1,"pv_univariate"], 1)
	checkEquals(res2[1,"pv_univariate"], 1)
	# merged has 2 phenotypes, but the permutations are intra-groups==> pv ==1
	checkEquals(res.merged[1,"pv_univariate"], 1)

	appended <- merged
	# artificially suppress the groups
	appended@phdata$gws <- NULL

	# appended has no more groups anm 2 phenotypes ==> pv < 1
	res.appended <-processGws(appended, bins, params)
	checkTrue(res.appended[1,"pv_univariate"] < 1)

}

test.processFiles <- function() {
	params <- parameters(
			types = c("univariate", "divariate"),
			seed = 0,
			verbosity = 0,
			nb_permutations = 100)

	res <- processFiles("extdata/ms1", params=params)

	b1 <- res[1,]
	checkEquals(b1$nb_snps, 1	)

	checkEqualsNumeric(b1$pv_univariate[1], 0.792079, tolerance=1e-6)
	checkEquals(b1$pv_divariate[1], 1)
	checkEqualsNumeric(res$fdr_univariate[1], 0.932665, tolerance=1e-6)

}

test.processGws <- function() {
	data("ms1")

	params <- parameters(seed = 0, verbosity = 0, nb_permutations = 100)
	res <- processGws(ms1_gws, ms1_bins, params)

	checkEqualsNumeric(res$pv_univariate[1], 0.792079, tolerance=1e-6)
	checkEqualsNumeric(res$fdr_univariate[1], 0.932665, tolerance=1e-6)

	checkEquals(dim(res)[1], 272)

	# test the phenotypes parameter
	res2 <- processGws(ms1_gws, ms1_bins, params, phenotypes=phdata(ms1_gws)$pop)
	checkEqualsNumeric(res$score_univariate, res2$score_univariate)

	# use different phenotype
	res3 <- processGws(ms1_gws, ms1_bins, params, phenotypes=rep(0, nids(ms1_gws)))
	checkTrue( ! isTRUE(all.equal.numeric(res$score_univariate, res3$score_univariate)) )
}

# use GenABEL dataset srdta
test.processGws.on.srdta <- function() {
	data(srdta)

	binsfile <-  "extdata/ms1.bins"
	bins <- readBins(binsfile)

	gws <- asGws(srdta, assignPopNATo = 0)

	params <- parameters(
			seed = 0,
			verbosity = 0,
			nb_permutations = 100)

	res <- processGws(gws, bins, params)

	checkEquals(dim(res)[1], 63)

}

# use GenABEL dataset srdta
test.processGws.no.bins <- function() {
	data("ms1")

	params <- parameters(seed = 0, verbosity = 0, nb_permutations = 1000)
	res <- processGws(ms1_gws, params=params)
	checkEqualsNumeric(nrow(res), nsnps(ms1_gws))


}


test.RcppRunNNBCOnFiles.badtype <- function() {
	f <- function() {
		params <- parameters(types = c("univariate", "titi"))
		res <- RcppRunNNBCOnFiles("extdata/ms1", parameters = params)
	}

	checkException(f(), silent = TRUE)
}

test.processTable <- function() {
	params <- parameters(
			seed = 0,
			nb_permutations = 10000)

	table2x4 <- matrix(c(100, 0, 20, 1,  80, 10, 10, 0), nrow=2,ncol=4, byrow=TRUE)

	res <- processTable(table2x4, params)

	checkEqualsNumeric(res$pvalues["univariate"], 0.000799920008, tolerance=1e-10)
}

test.processTable_on_divariate_table <- function() {
	params <- parameters(
			seed = 0,
			nb_permutations = 10000)

#	params$types <- "univariate"

	table2x4x4 <- c(277,250,5,2,0,0,0,1,    107,106,3,3,0,0,0,0,   8, 12, 1, 1,0,0,0,0,   1,1,0,0,0,0,0,0)
	dim(table2x4x4) <- c(2,4,4)
#	print(table[1,,])
#	print(table[2,,])
	res <- processTable(table2x4x4, params)

	checkEqualsNumeric(res$pvalues["divariate"], 0.8421157884, tolerance=1e-10)

}
