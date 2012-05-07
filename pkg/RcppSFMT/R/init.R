.onLoad <- function(lib, pkg)
{
	library.dynam("rngSFMT2", pkg, lib)
}

.onUnload <- function(libpath) {
	library.dynam.unload("rngSFMT2", libpath)
}

