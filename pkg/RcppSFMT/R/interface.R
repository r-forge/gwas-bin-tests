MERSENNE_EXPONENT <- function() {
	.Call("MERSENNE_EXPONENT",PACKAGE = "RcppSFMT")
}

SFMT_ID <- function() {
	.Call("SFMT_ID",PACKAGE = "RcppSFMT")
}

rand32 <- function(nb=1, seeds=1) {
	.Call("W_fill_array32_seeds", nb, seeds, PACKAGE = "RcppSFMT")	
}

