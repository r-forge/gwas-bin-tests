library(rngVerify)
n <- 50 # number of parallel streams
m <- 10000 # total number of random numbers generated
# sequence of m indices of streams from {0, ..., n-1}, small indices are more frequent
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
x <- testInterleaving(n, as.integer(ind))
outInter <- unlist(lapply(split(x, ind), sumDouble))

.Random.seed <- keep
outSeq <- testSequential(n, len)

.Random.seed <- keep
outThreaded <- testThreaded(n, len)

err1 <- max(abs(outSeq - outInter))
err2 <- max(abs(outSeq - outThreaded))

err <- c(err1, err2)
print(err)

# The following test requires compiler option "-ffloat-store" on
# some 32-bit Intel architectures, see rngVerify/src/Makevars.
stopifnot(err == 0)

