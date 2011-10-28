#' dynamic loading
#' @param lib
#' @param pkg
#' @import GenABEL
#' @author	karl

.onLoad <- function(lib, pkg) {
	library.dynam("GWASBinTests", pkg, lib )
}

.onUnload <- function (libpath) {
	cat(".onUnload, libpath=", libpath)
	library.dynam.unload("GWASBinTests", libpath)
}