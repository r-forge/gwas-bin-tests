MRG32k5a <- function(state, n)
{
    norm <- 2.3283163396834613e-10
    m1 <- as.bigz(2^32 - 18269)
    m2 <- as.bigz(2^32 - 32969)
    a12 <- as.bigz( 1154721)
    a14 <- as.bigz( 1739991)
    a15 <- as.bigz(-1108499)
    a21 <- as.bigz( 1776413)
    a23 <- as.bigz(  865203)
    a25 <- as.bigz(-1641052)
    x <- rep(NA, times=n)
    state <- as.bigz(state)
    for (i in seq.int(length.out=n)) {
        p1 <- mod.bigz(a12*state[4]  + a14*state[2] + a15*state[1], m1)
        p2 <- mod.bigz(a21*state[10] + a23*state[8] + a25*state[6], m2)
        state <- c(state[2:5], p1, state[7:10], p2)
        x[i] <- as.double(mod.bigz(p1 - p2, m1)) * norm
    }
    x
}

