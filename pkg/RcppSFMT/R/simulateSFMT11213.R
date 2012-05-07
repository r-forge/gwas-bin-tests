simulateSFMT11213 <- function(state, n)
{
    # 128-bit integers are represented by vectors of bits, where the most significant
    # bit is at index 1 and the least significant bit is at index 128

    hexa <- function(s)
    {
        dig <- c(0:9, letters[1:6])
        state <- match(strsplit(s, "")[[1]], dig) - 1
        c(sapply(state, FUN = function(y) { floor(y / 2^(3:0)) %% 2 }))
    }

    A <- function(w)
    {
        s <- rep(0, times=128)
        s[1:104] <- w[25:128]
        w + s
    }

    BShift <- c(1:25, 1:25 + 32, 1:25 + 64, 1:25 + 96)
    BMask <- hexa(paste("7fffdbfd", "dfdfbfff", "ffffffef", "effff7fb", sep=""))
    B <- function(w)
    {
        s <- rep(0, times=128)
        s[BShift + 7] <- w[BShift]
        s * BMask
    }

    C <- function(w)
    {
        s <- rep(0, times=128)
        s[25:128] <- w[1:104]
        s
    }

    DShift <- c(15:32, 15:32 + 32, 15:32 + 64, 15:32 + 96)
    D <- function(w)
    {
        s <- rep(0, times=128)
        s[DShift - 14] <- w[DShift]
        s
    }

    recursion <- function(S)
    {
        X <- A(S[1, ]) + B(S[69, ]) + C(S[87, ]) + D(S[88, ])
        S[1:87, ] <- S[2:88, ]
        S[88, ] <- X %% 2
        S
    }

    stopifnot(length(state) == 352)
    X <- sapply(state, FUN = function(y) { floor(y / 2^(0:31)) %% 2 } )
    S <- t(matrix(X, nrow=nrow(X)*4, ncol=ncol(X)/4))[, 128:1]
    nn <- ceiling(n/4)
    y <- rep(NA, times=4*nn)
    for (i in seq.int(length.out=nn)) {
        S <- recursion(S)
        a <- matrix(S[88, ], nrow=32, ncol=4)
        y[4*(i-1) + 1:4] <- rev(rbind(2^(31:0)) %*% a)
    }
    y[1:n]
}

testSFMT11213 <- function(x)
{
    stopifnot(length(x) > 352)
    y <- simulateSFMT11213(x[1:352], length(x) - 352)
    table(y == x[-(1:352)], useNA="ifany")
}

