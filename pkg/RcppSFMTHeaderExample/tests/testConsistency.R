library(rngVerify)
n <- 50 # number of parallel streams
m <- 10000 # total number of random numbers generated
# sequence of m indices of streams from 1:n, small indices are more frequent
ind <- factor(ceiling(n*runif(m)^2), levels=1:n)
len <- c(table(ind))

# replace sum(), which works in long double, by sumDouble
# to make the error smaller
sumDouble <- function(xx)
{
	s <- 0
	for (x in xx) s <- s + x
	s
}

runif(1)

keep <- .Random.seed
interSeq <- tryInterleavingSeq(as.integer(ind))

.Random.seed <- keep
interThr <- tryInterleavingThr(as.integer(ind))

stopifnot(identical(interSeq, interThr))

.Random.seed <- keep
sumSeq <- trySumSeq(len)

.Random.seed <- keep
sumThr <- trySumThr(len)

err1 <- max(abs(sumSeq - sumThr))
print(err1)
stopifnot(err1 < 1e-10)

sumThrSplit <- unlist(lapply(split(interThr, ind), sumDouble))

err2 <- max(abs(sumThrSplit - sumThr))
print(err2)
stopifnot(err2 < 1e-10)

