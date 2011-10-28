#library(RUnit)

test.fetchGenotypes <- function(parameter) {
	data(srdta)
	g <- fetchGenotypes(srdta,1)
	checkTrue(length(g) == 2500 && sum(g == 0) == 1792)

	# index range
	index <- 10:20
	g <- fetchGenotypes(srdta, index)
	checkTrue( all(dim(g) == c(2500,length(index))) )

	# by names
	names <- snp.names(srdta)[index]
	g2 <- fetchGenotypes(srdta, index)
	checkTrue( all(g==g2) )

	# all
	g <- fetchGenotypes(srdta)
	checkTrue( all(dim(g) == c(2500,nsnps(srdta))) )

}

test.fetchGenotypesAsList <- function() {
	data(srdta)
	g <- fetchGenotypesAsList(srdta,1)
	checkTrue(length(g) == 1 && length(g[[1]]) == 2500 && sum(g[[1]] == 0) == 1792)

	# index range
	index <- 10:20
	g <- fetchGenotypesAsList(srdta, index)
	checkTrue( length(g) == length(index)  && length(g$rs150) == 2500)

	# by names
	names <- snp.names(srdta)[index]
	g2 <- fetchGenotypesAsList(srdta, index)
	checkEquals(g, g2)

	# all
	g <- fetchGenotypesAsList(srdta)
	checkEqualsNumeric( length(g), nsnps(srdta) )

}

test.missingRateFromGenotypes <- function() {
	data(srdta)
	beta <- missingRateFromGenotypes( fetchGenotypes(srdta, 1) )
	checkEqualsNumeric(beta, 116/2500)
}

test.convertGenotypes <- function(parameter) {
	# check type
	checkException( convertGenotypes(c("coucou", " hello")), silent = TRUE )

	# check values
	checkException( convertGenotypes(1:10), silent = TRUE )

	# check empty vector
	checkException( convertGenotypes(integer(0)), silent = TRUE )

	data(srdta)
	gs <- convertGenotypes( fetchGenotypes(srdta, 1:20) )
	checkTrue( all( dim(gs) == c(2500,20)) )

}