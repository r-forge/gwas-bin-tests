# library(RUnit)
# library(GWASBinTests, lib.loc = "..truetest/")

test.gwsForBin <- function() {
	data(srdta)

	# empty gws
	gws <- gwsForBin(srdta, 30,1,1)
	checkTrue( is.null(gws) )

	start <- map(srdta[,10])
	end <- map(srdta[,20])

	gws <- gwsForBin(srdta, 1, start, end)
	checkTrue( nsnps(gws) == 11 && nids(gws) == nids(srdta))

}

test.makeBinsForSnps <- function() {
	data(srdta)

	bins <- makeBinsForSnps(srdta)

	checkTrue(nrow(bins) == nsnps(srdta))
	chrs <- convertChromosomeNames(bins$chr_name)
	checkTrue( ! is.unsorted(chrs) )

	by_bin <- apply(bins[40:49,],1,function(b) {
		gwsForBin(srdta, b["chr_name"], b["start"], b["end"])
	})
	nsnp_by_bin <- sapply(by_bin, nsnps)
	checkTrue(  all( nsnp_by_bin == 1))

}

