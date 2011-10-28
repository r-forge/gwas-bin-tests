#library(RUnit)

test.simpleGenotypeErrorModel <- function() {
	alpha <- 0.05
	beta <- 0.01
	m <- simpleGenotypeErrorModel(alpha, beta)

	checkTrue( all.equal.numeric( diag(m), rep( (1-beta)*(1-alpha), 3), check.names=FALSE, check.attributes = FALSE ) )
	checkEqualsNumeric( m["aa","Aa"], (1-beta)*alpha/2 )
	checkEqualsNumeric(m["ZZ",], rep(beta,3) )

	m <- simpleGenotypeErrorModel(0, 0)
	checkTrue( all.equal.numeric( diag(m), rep( 1, 3),, check.names=FALSE, check.attributes = FALSE ) )
	diag(m) <- 0
	checkTrue( all(m ==0 | is.na(m)) )

}

