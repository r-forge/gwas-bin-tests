
#library(RUnit)

test.saveGws <- function() {
	data(srdta)
	gws <- asGws(srdta[1:100,], assignPopNATo = 0)
	file <- tempfile("gws")

	saveGws(gws, file)

	new <- readGws(file)

	checkTrue( new@gtdata@nsnps == gws@gtdata@nsnps)
	checkEquals(gws@phdata$ids, new@phdata$ids)

}
