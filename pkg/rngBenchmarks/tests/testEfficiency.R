library(rngVerify)
n <- 10 # number of parallel streams
m <- 5e7 # number of random numbers generated from each stream

runif(1)
keep <- .Random.seed
tThreaded <- system.time(outThreaded <- testThreaded(n, rep(m, times=n)))

.Random.seed <- keep
tSeq <- system.time(outSeq <- testSequential(n, rep(m, times=n)))

err <- abs(outSeq - outThreaded)
print(err)
stopifnot(err == 0)

print(rbind(tSeq, tThreaded))

tSeq[3]/tThreaded[3]