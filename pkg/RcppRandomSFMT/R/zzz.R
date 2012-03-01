#' dynamic loading
#' @param lib
#' @param pkg
#' @import RcppRandomSFMT
#' @author	karl

.onLoad <- function(lib, pkg) {
	library.dynam("RcppRandomSFMT", pkg, lib )
}

.onUnload <- function (libpath) {
	cat(".onUnload, libpath=", libpath)
	library.dynam.unload("RcppRandomSFMT", libpath)
}