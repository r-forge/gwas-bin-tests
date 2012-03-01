library(rngBenchmarks)
library(car)
library(plyr)
library(lattice)

types <- c("pureR_sequential_vectorized", "random_r", "lecuyer", "rng_openmp", "sfmt")
# check reproducibility
n <- 1e6
res <- sapply(types, function(type) {
			cat("type=", type, "\n")
			pi1 <- pi_rng_threaded(n, threads=1, seed=2, type=type)
			cat
			pi2 <- pi_rng_threaded(n, threads=1, seed=2, type=type)
			abs(pi1-pi2)
})
stopifnot(all( res == 0 ))

fast_types <- c("random_r", "lecuyer", "rng_openmp", "sfmt")
#fast_types <- c("rng_openmp", "sfmt")
nbs <- c(1e7, 1e8, 1e9)
#nbs <- c(1e5, 1e6)
nb_threads <- c(1,2,4,8)

dat <- expand.grid( list(type=fast_types, threads=nb_threads,n=nbs))
dat$type <- as.character(dat$type)

# subset(dat, type=="sfmt" & threads==1 & n==1e05)
res <- ddply( dat, .variables = colnames(dat), function(tab) {

			tim <- system.time(pih <- pi_rng_threaded(n=tab$n, threads=tab$threads, seed=3, type=tab$type), gcFirst = TRUE)
			tab$error <- pi_error(pih)
			tab$time <- tim[3]
			tab

		})
some(res)
res$time_per_chunk <- res$time / res$n * 10000;
res$time_per_chunk_per_thread <- res$time / res$n *res$threads * 10000;

xyplot( time ~ n| as.factor(threads), type="b" ,groups=type, auto.key=TRUE,res)

xyplot( time_per_chunk ~ threads, type="b" ,groups=type, auto.key=TRUE,res)
xyplot( time_per_chunk_per_thread ~ as.factor(type), groups=threads, type="b", auto.key=TRUE,res)

xyplot( error ~ n| as.factor(n), type="b" ,groups=type, auto.key=TRUE,res)


# more benchs for 8 threads
fast_types <- c("random_r", "lecuyer", "rng_openmp", "sfmt")
nbs <- c(1e7, 1e8, 1e9, 1e10)
nb_threads <- 8
dat <- expand.grid( list(type=fast_types, threads=nb_threads,n=nbs))
dat$type <- as.character(dat$type)
res8 <- ddply( dat, .variables = colnames(dat), function(tab) {

			tim <- system.time(pih <- pi_rng_threaded(n=tab$n, threads=tab$threads, seed=3, type=tab$type), gcFirst = TRUE)
			tab$error <- pi_error(pih)
			tab$time <- tim[3]
			tab
		})


xyplot( time ~ n| as.factor(threads), type="b" ,groups=type, auto.key=TRUE,res8)
xyplot( error ~ n| as.factor(n), type="b" ,groups=type, auto.key=TRUE,res)

# save results
previous_results <- res


#    user  system elapsed
#   0.808   0.000   0.810


pi_error(pirl)
# [1] 9.273844e-05


setwd("~/workspace/rngBenchmarks")
library(rngOpenMP, lib.loc="~/R")
library(RcppRandomSFMT, lib.loc="~/R")
library(rngBenchmarks, lib.loc="..truetest/")

unloadNamespace("rngBenchmarks")