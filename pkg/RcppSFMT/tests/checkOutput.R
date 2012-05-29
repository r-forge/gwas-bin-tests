#library(RcppSFMT, lib.loc="../tmpsse/")
library(RcppSFMT)

mexp <- new(SFMT)$get_mersenne_exponent()
cat("SFMT MERSENNE EXPONENT =", mexp, "\n")
file <- paste("SFMT.", mexp, ".out.txt", sep="")
path <- system.file("include", "SFMT-src", file, package="RcppSFMT")

if ( path == "" )
	warning("Could not find file ", file)
stopifnot(path != "")

parseLines <- function(lines) {
	v <- gsub('\\s+$', '', lines, perl = TRUE) # trailing spaces
	v <- gsub('^\\s+', '', v, perl = TRUE) # front spaces
	v <- as.numeric(unlist(strsplit(v, '\\s+', perl=TRUE)))
	v
}

con <- file(path, "r")

dummy <- readLines(con, n=3)  # skip header
# read 1000 values, 5 by line => 200 lines

init_gen_rand <- parseLines( readLines(con, n=200) )
stopifnot(length(init_gen_rand) == 1000 )

dummy <- readLines(con, n=2)  # skip header2
init_by_array <- parseLines( readLines(con, n=200) )
stopifnot(length(init_by_array) == 1000 )

close(con)

# generate random numbers
rng <- new(SFMT, 1234)
rands <- rng$unsignedRandomIntegers(1000)
# check them
stopifnot( all.equal(rands, init_gen_rand) )

# use init_by_array, i.e. multiple seeds
ini <-  c(0x1234, 0x5678, 0x9abc, 0xdef0)
rng$set_seed_array(ini)
rands2 <- rng$unsignedRandomIntegers( 1000 )
stopifnot( all.equal(rands2, init_by_array) )