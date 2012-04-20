simulateSFMT19937 <- function(state, n)
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
        s[1:120] <- w[9:128]
        w + s
    }

    BShift <- c(1:21, 1:21 + 32, 1:21 + 64, 1:21 + 96)
    BMask <- hexa(paste("bffffff6", "bffaffff", "ddfecb7f", "dfffffef", sep=""))
    B <- function(w)
    {
        s <- rep(0, times=128)
        s[BShift + 11] <- w[BShift]
        s * BMask
    }

    C <- function(w)
    {
        s <- rep(0, times=128)
        s[9:128] <- w[1:120]
        s
    }

    DShift <- c(19:32, 19:32 + 32, 19:32 + 64, 19:32 + 96)
    D <- function(w)
    {
        s <- rep(0, times=128)
        s[DShift - 18] <- w[DShift]
        s
    }

    recursion <- function(S)
    {
        X <- A(S[1, ]) + B(S[123, ]) + C(S[155, ]) + D(S[156, ])
        S[1:155, ] <- S[2:156, ]
        S[156, ] <- X %% 2
        S
    }

    stopifnot(length(state) == 624)
    X <- sapply(state, FUN = function(y) { floor(y / 2^(0:31)) %% 2 } )
    S <- t(matrix(X, nrow=nrow(X)*4, ncol=ncol(X)/4))[, 128:1]
    nn <- ceiling(n/4)
    y <- rep(NA, times=4*nn)
    for (i in seq.int(length.out=nn)) {
        S <- recursion(S)
        a <- matrix(S[156, ], nrow=32, ncol=4)
        y[4*(i-1) + 1:4] <- rev(rbind(2^(31:0)) %*% a)
    }
    y[1:n]
}

testSFMT19937 <- function(x)
{
    stopifnot(length(x) > 624)
    y <- simulateSFMT19937(x[1:624], length(x) - 624)
    table(y == x[-(1:624)], useNA="ifany")
}

