testInterleaving <- function(n, ind)
{
    m <- length(ind)
    out <- .C("testInterleaving",
        as.integer(n),
        as.integer(m),
        as.integer(ind - 1), # convert to C indexing
        out=double(m),
        PACKAGE="rngBenchmarks")
    out$out
}

testSequential <- function(n, len)
{
    out <- .C("testSequential",
        as.integer(n),
        as.integer(len),
        out=double(n),
        PACKAGE="rngBenchmarks")
    out$out
}

testThreaded <- function(n, len)
{
    out <- .C("testThreaded",
        as.integer(n),
        as.integer(len),
        out=double(n),
        PACKAGE="rngBenchmarks")
    out$out
}

