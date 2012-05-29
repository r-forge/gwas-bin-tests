
.onLoad <- function(lib, pkg)
{
	require("methods", character=TRUE, quietly=TRUE)
	loadRcppModules()
#	library.dynam("RcppSFMT", pkg, lib)
}

.onUnload <- function(libpath) {
	library.dynam.unload("RcppSFMT", libpath)
}

