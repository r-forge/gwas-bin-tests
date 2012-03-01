#' dynamic loading
#' @param lib
#' @param pkg
#' @import rngBenchmarks
#' @author	karl

.onLoad <- function(lib, pkg) {
	library.dynam("rngBenchmarks", pkg, lib )
}

.onUnload <- function (libpath) {
	cat(".onUnload, libpath=", libpath)
	library.dynam.unload("rngBenchmarks", libpath)
}