trySetRngParams <- function(initType)
{
    .C("trySetRngParams",
        as.integer(initType),
        PACKAGE="rngVerify")
    invisible(NULL)
}

tryGetStreams <- function(n, m)
{
    out <- .C("tryGetStreams",
        as.integer(n),
        as.integer(m),
        out=double(n*m),
        PACKAGE="rngVerify")$out
    matrix(out, nrow=n, ncol=m)
}

tryInterleavingSeq <- function(ind)
{
    m <- length(ind)
    out <- .C("tryInterleavingSeq",
        as.integer(max(ind)),
        as.integer(m),
        as.integer(ind - 1), # convert to C indexing
        out=double(m),
        PACKAGE="rngVerify")
    out$out
}

tryInterleavingThr <- function(ind, nbThreads=8)
{
    m <- length(ind)
    out <- .C("tryInterleavingThr",
        as.integer(max(ind)),
        as.integer(m),
        as.integer(ind - 1), # convert to C indexing
        out=double(m),
        as.integer(nbThreads),
        PACKAGE="rngVerify")
    out$out
}

trySumSeq <- function(len)
{
    n <- length(len)
    out <- .C("trySumSequential",
        as.integer(n),
        as.integer(len),
        out=double(n),
        PACKAGE="rngVerify")
    out$out
}

trySumThr <- function(len, nbThreads=8)
{
    n <- length(len)
    out <- .C("trySumThreaded",
        as.integer(n),
        as.integer(len),
        out=double(n),
        as.integer(nbThreads),
        PACKAGE="rngVerify")
    out$out
}

testPi <- function(n, nbThreads=8)
{
    stopifnot(n >= 10000)
    out <- .C("pi_rng_threaded",
        as.integer(n),
        out=double(1),
        as.integer(nbThreads),
        PACKAGE="rngVerify")
    out$out
}

