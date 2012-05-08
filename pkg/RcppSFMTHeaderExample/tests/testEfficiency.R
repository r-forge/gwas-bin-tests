library(rngVerify)
n <- 10 # number of parallel streams
m <- 2e6 # number of random numbers generated from each stream

runif(1)
keep <- .Random.seed
timeThr <- system.time(outThr <- trySumThr(rep(m, times=n)))

.Random.seed <- keep
timeSeq <- system.time(outSeq <- trySumSeq(rep(m, times=n)))

err <- abs(outSeq - outThr)
print(err)
stopifnot(err < 1e-10)

print(rbind(timeSeq, timeThr))

