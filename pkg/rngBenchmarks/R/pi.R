# Monte-Carlo PI computation implementations
#
###############################################################################

rng_openmp <- function(n , threads, seed){
	set.seed(seed);
	.Call("_pi_rngOpenMP_threaded",  n , threads, seed, PACKAGE="rngBenchmarks")
}

pi_rng_threaded <- function(n, threads, type, seed = 1) {
	switch(type,
		sfmt = .Call("_pi_rng_sfmt_threaded",  n , threads, seed, PACKAGE="rngBenchmarks")
		,lecuyer=.Call("_pi_rng_lecuyer_threaded",  n , threads, seed, PACKAGE="rngBenchmarks")
		,random_r=.Call("_pi_rng_random_r_threaded",  n , threads, seed, PACKAGE="rngBenchmarks")
		,pureR_sequential = rngBenchmarks:::pi_pureR_sequential(n, seed)
		,pureR_sequential_vectorized = rngBenchmarks:::pi_pureR_sequential_vectorized(n, seed)
		,rng_openmp=rng_openmp(n, threads, seed)
	)
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