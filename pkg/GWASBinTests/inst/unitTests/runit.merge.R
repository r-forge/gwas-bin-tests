#library(GenABEL, lib.loc = "~/R")
#library(RUnit)
#library(GWASBinTests, lib.loc = "..truetest/")

test.mergeGws <- function() {
	data(srdta)

	# make srdta compatible with GWASBinTests
	srdta@phdata$pop <- srdta@phdata$bt
	gws <- asGws(srdta, assignPopNATo = 0)

	# distinct samples

	i1 <- 1:4
	i2 <- 10:12
	x1 <- gws[i1,]
	x2 <- gws[i2,]
	x <-mergeGws(x1, x2)
	GWASBinTests:::.checkPhenoData( x@phdata )

	xi1 <- 1:length(i1)
	xi2 <- 1:length(i2) + length(xi1)

	for ( var in c("id", "sex", "pop") ) {
		checkEquals(x1@phdata[[var]], x@phdata[[var]][xi1])
		checkEquals(x2@phdata[[var]], x@phdata[[var]][xi2])
	}

	xgws <- x@phdata[["gws"]]
	checkTrue( all(xgws[xi1] == 1) )
	checkTrue( all(xgws[xi2] == 2) )

	checkTrue(is.integer(x@phdata$sex))
	checkTrue(is.integer(x@phdata$pop))
	checkTrue(is.integer(x@phdata$gws))

	# merge with common ids should fail
	checkException(mergeGws(x1, x1), silent=TRUE )


	i3 <- 20:30
	x3 <- gws[i3,]
	xx <- mergeGws(x, x3)
	checkEquals(sort(unique(xx@phdata$gws)), 1:3)

	# merge on same sample with distinct ids
	x1 <- gws[2,]
	x2 <- gws[2,]
	# change id
	x1@gtdata@idnames <- "fromx1"
	x1@phdata$id <- x1@gtdata@idnames
	x <- mergeGws(x1, x2)


	checkTrue( nsnps(x) == 833 && nids(x) == 2 )

	## merge on same samples, different SNPs
	## not allowed for now
	x1 <- gws[1:4,1]
	x2 <- gws[5:10, 2]
	# no SNP in common
	checkException( mergeGws(x1, x2, intersected_snps_only = TRUE), silent = TRUE )

	merge.snp.data(x1@gtdata, x2@gtdata, intersected_snps_only = FALSE )
	x <- mergeGws(x1, x2, intersected_snps_only = FALSE)
	# genotypes for SNP2 replaced by NA
	checkTrue( all(is.na(as.numeric(x[1:4, 2]))) )
	checkEqualsNumeric(as.numeric(x[1:4, 1]), as.numeric(x1) )
	# idem for SNP1
	checkTrue( all(is.na(as.numeric(x[5:10, 1]))) )
	checkEqualsNumeric(as.numeric(x[5:10, 2]), as.numeric(x2) )

}

