test_SFMT <- function(seed, nb=1) {
	.C("sfmt_test", as.integer(seed),as.integer(nb),
					out=double(1), PACKAGE = "RcppSFMT")$out
}


SFMT_printID <- function() {
	.Call("sfmt_printid",PACKAGE = "RcppSFMT")
}

initSFMT <- function(seed=NULL, vseed=NULL)
{
    useArray <- is.null(seed)
    if (useArray) {
        seed <- vseed
    } else {
        stopifnot(length(seed) == 1)
    }
	stopifnot(seed >= 0)
	stopifnot(seed < 2^32)
	stopifnot(floor(seed) == seed)
    seed[seed >= 2^31] <- seed[seed >= 2^31] - 2^32
    seed[seed == -2^31] <- NA

    if (useArray) {
		.C("sfmt_init_array",
			length(seed),
			as.integer(seed),
			PACKAGE="rngSFMT2")
    } else {
		.C("sfmt_init_single",
			as.integer(seed),
			PACKAGE="rngSFMT2")
    }
    invisible(NULL)
}

generIntSFMT <- function(n)
{
    .C("sfmt_gener_int",
        as.integer(n),
        out=double(n),
        PACKAGE="rngSFMT2")$out
}

generRealSFMT <- function(n)
{
    .C("sfmt_gener_real",
        as.integer(n),
        out=double(n),
        PACKAGE="rngSFMT2")$out
}

deleteSFMT <- function()
{
    .C("sfmt_delete",
        PACKAGE="rngSFMT2")
    invisible(NULL)
}

formatOutput <- function(useArray)
{
    if (useArray) {
        cat("init_by_array\n")
        initSFMT(vseed=as.numeric(c("0x1234", "0x5678", "0x9abc", "0xdef0")))
        x <- generIntSFMT(1000)
    } else {
        cat("init_gen_rand\n")
        initSFMT(seed=1234)
        x <- generIntSFMT(1000)
    }
    for (i in 1:200) {
        for (j in 1:5) {
            cat(sprintf("%10.0f ", x[5*(i-1) + j]))
        }
        cat("\n")
    }
    invisible(NULL)
}

