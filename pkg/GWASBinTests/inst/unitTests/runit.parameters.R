
#library(RUnit)

test.default.parameters <- function() {
	params <- parameters()
	checkTrue( params["max_error_rate"] < 1)
}

#test.parameters.rcpp  <- function() {
#	params <- parameters(nb_permutations=7, verbosity=1,
#				use_affymetrix_model=1, confidence=0.12, threads=6, types = c("univariate", "divariate"))
#	.Call("RcppParamsExample", params)
#
#}

