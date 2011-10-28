
#library(RUnit)

.checkLivesOk <- function(expr) {
	error <- FALSE
	tryCatch(expr, error = function(e) { error <<- TRUE} )
	checkTrue( !error )
}

test.readGws <- function() {
	file <- "extdata/ms1"
	gws <- readGws(file, verbose=FALSE)
	checkTrue( gws@gtdata@nsnps == 1000)
}

test.readBins <- function() {
	file <- "extdata/ms1.bins"
	bins <- readBins(file)
	checkTrue( is.data.frame(bins) && length(bins[,1]) == 21549 )

	# check it it sorted
	chrs <- convertChromosomeNames(bins$chr_name)
	checkTrue( ! is.unsorted(chrs) )
}

test.asGws <- function() {
	require(GenABEL)
	data(srdta)

	# this should fail because of the NAs in phdata$bt
	checkException( gws <- asGws( srdta ), silent = TRUE)

	# this should now work
	.checkLivesOk( gws <- asGws( srdta, assignPopNATo=0) )

	.checkLivesOk(GWASBinTests:::.checkPhenoData( gws@phdata ))

}

test.readPlinkTransposedData <- function() {
	file <- "extdata/plink_small"
	gws <- readPlinkTransposedData(file, verbose = FALSE)
	#print(gws)
}