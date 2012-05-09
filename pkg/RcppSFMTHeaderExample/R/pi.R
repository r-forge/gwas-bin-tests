# Monte-Carlo PI computation implementations
#
###############################################################################


pi_rng_threaded <- function(n, threads, seed = 1, chunk_size=10000) {
	.Call("_pi_rng_threaded",  n , threads, seed, chunk_size,PACKAGE="RcppSFMTHeaderExample")
}

pi_pureR_sequential <- function(n, seed) {
	set.seed(seed)
	inside <- 0
	for (i in 1:n) {
		p <- runif(2, min=-1, max=1)
		if ( sum(p*p)<= 1)
			inside <- inside + 1

	}
	pi_hat <- inside/n*4
}

pi_pureR_sequential_vectorized <- function(n, seed) {
	set.seed(seed)
	inside <- 0
	v <- runif(2*n, min=-1, max=1)
	v <- v*v
	p <- v[1:(2*n-1)] + v[2:(2*n)]
	for (i in 1:n) {
		if ( p[2*i-1] <= 1)
			inside <- inside + 1
	}
	pi_hat <- inside/n*4
}


pi_error <- function(pih) {
	abs(pih-pi)/pi
}