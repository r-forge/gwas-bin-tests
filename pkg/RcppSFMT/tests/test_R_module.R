library(RcppSFMT)

# ===== constructors ======

# default constructor => used R seed
set.seed(1)
r <- new(SFMT)

n1 <- r$nextRandomInteger()
set.seed(1)
r$set_seeds_from_R()
n2 <- r$nextRandomInteger()
stopifnot( n1 == n2)

# constructor with seed
r <- new(SFMT, 1234)
n1 <- r$nextRandomInteger()
r$set_seed(1234)
n2 <- r$nextRandomInteger()
stopifnot( n1 == n2)

# constructor with 2 seeds
r <- new(SFMT, 1234, 3456)
n1 <- r$nextRandomInteger()
r$set_seeds(1234, 3456)
n2 <- r$nextRandomInteger()
stopifnot( n1 == n2)

# id and mersenne exponent
r <- new(SFMT)
mexp <- r$get_mersenne_exponent()
stopifnot( mexp > 0 )

id <-  r$id() 
stopifnot( nchar(id) > 0)

# ===== set_seed_array and set_seeds======
r <- new(SFMT)
r$set_seeds(1234, 3456)
n1 <- r$nextRandomInteger()
r$set_seed_array( c(1234, 3456))
n2 <- r$nextRandomInteger()
stopifnot( n1 == n2)

# randomIntegers
r <- new(SFMT, 1234)
ris1 <-  sapply(1:397, function(x) { r$nextRandomInteger( ) })
r$set_seed(1234)
ris2 <-  r$randomIntegers( 397 )
stopifnot(  all.equal(ris1, ris2) )

# unsignedRandomIntegers
r <- new(SFMT, 1234)
ris1 <-  sapply(1:397, function(x) { r$nextUnsignedRandomInteger( ) })
r$set_seed(1234)
ris2 <-  r$unsignedRandomIntegers( 397 )
stopifnot(  all.equal(ris1, ris2) )

# randomDoubles
r <- new(SFMT, 1234)
ris1 <-  sapply(1:397, function(x) { r$nextRandomDouble( ) })
r$set_seed(1234)
ris2 <-  r$randomDoubles( 397 )
stopifnot(  all.equal(ris1, ris2) )

stopifnot( all( ris1 > 0 & ris1 < 1 )  )




