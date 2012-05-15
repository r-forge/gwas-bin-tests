MERSENNE_EXPONENT <- function() {
	.Call("MERSENNE_EXPONENT",PACKAGE = "RcppSFMT")
}

rand32 <- function(nb=1, seed=1) {
	.Call("W_fill_array32", nb, seed, PACKAGE = "RcppSFMT")	
}