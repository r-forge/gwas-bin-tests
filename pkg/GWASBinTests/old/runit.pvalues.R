#library(RUnit)

test.permuteByGroup <- function() {

	# test permutation keeps group properties
	values <- cbind(1:20)
	indices_list <- list(g2=which(values > 5 & values%%2==1), g1=which(values%%2==0))
	for(i in 1:10) {
		values <- permuteByGroup(values, indices_list)
		checkEquals(which(values > 5 & values%%2==1), indices_list$g2)
	}

	# test one group
	values <- cbind(1:10)
	indices_list <- list(1:2)
	values <- permuteByGroup(values, indices_list)
	checkEquals(values[3:10], 3:10)

	# test one singleton
	values <- cbind(1:10)
	indices_list <- list(1)
	values <- permuteByGroup(values, indices_list)
	checkEquals(values, cbind(1:10))

	# test empty group
	values <- cbind(1:10)
	indices_list <- list(g=c())
	values <- permuteByGroup(values, indices_list)
	checkEquals(values, cbind(1:10))
}

test.permuteVariables <- function() {

	values <- 1:10
	df <- data.frame(v=values,g1=as.integer(values > 5 & values%%2==1), g2=as.integer(values%%2==0))
	res <- sapply(1:10, function(t) {
		x <- permuteVariables(df)
		min(x[df$g1==1,1])
	})
	checkTrue( all(res==7))

	## test: no-covariable
	values <- 1:100
	df <- data.frame(v=values)
	res <- sapply(1:100, function(t) {
				x <- permuteVariables(df)
				x[1,1]
	})
	m <- mean(res)
	checkTrue( abs(m-50) < 10)
}

#test.estimatePvalue <- function() {
#	values <- 1:10
#	df <- data.frame( S=c(rep(1,5), rep(0,5)))
#
#
#	# this score is minimial for the observed configuration
#	score <- function(vars) {
#		sum( values[vars$S==1])*2 + sum(values[vars$S==0])
#	}
#	res <- estimatePvalue(df, 100, score)
#	checkEquals(res$pvalue, 1)
#
#	# this score is maximal for the observed configuration
#	score2 <- function(vars) {
#		sum( values[vars$S==1]) + sum(values[vars$S==0])*2
#	}
#	res <- estimatePvalue(df, 100, score2)
#	checkEquals(res$pvalue, 1/101)
#
#	# this score is pefectly random
#	score3 <- function(vars) {
#		vars$S[1]
#	}
#	res <- estimatePvalue(df, 100, score3)
#	checkTrue(res$pvalue > 0.1 && res$pvalue < 0.9)
#
#}
