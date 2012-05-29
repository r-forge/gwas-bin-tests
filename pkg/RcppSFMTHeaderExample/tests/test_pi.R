#library(RcppSFMT, lib.loc="tmpsse/")
#library(RcppSFMTHeaderExample, lib.loc="tmpsse/")
library(RcppSFMTHeaderExample)

cat('ID=',new(SFMT)$id(), "\n")

cat("using one thread\n")
system.time( pi1 <- RcppSFMTHeaderExample:::pi_rng_threaded(n=1e8, threads=1, seed=1) )

cat("using 4 threads\n")
system.time( pi2 <- RcppSFMTHeaderExample:::pi_rng_threaded(n=1e8, threads=4, seed=1) )

stopifnot(pi1 == pi2)

stopifnot( RcppSFMTHeaderExample:::pi_error(pi1) < 1e-3 )


