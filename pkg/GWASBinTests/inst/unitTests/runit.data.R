#library(RUnit)

test.ms1 <- function() {
	files <- paste("extdata/ms1.",c("gag", "gap", "bins"), sep="")
	checkTrue( all( file.exists(files) ) )

	rm( list = ls() )
	data("ms1")
	checkTrue( ! is.null(ms1_gws) && nsnps(ms1_gws) == 1000 && nids(ms1_gws) == 778 )
	checkTrue( ! is.null(ms1_bins) && dim(ms1_bins) == c(21549, 4) )

}

test.plink_small <- function() {
	files <- paste("extdata/plink_small.",c("tfam", "tped"), sep="")
	checkTrue( all(file.exists(files) )  )
}