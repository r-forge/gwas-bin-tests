## This file has been shameless copied and adapted from RcppGSL/R/inline.R file
## See below for the original license and copyright information
##
## Copyright (C)       2010 Dirk Eddelbuettel and Romain Francois
##
## This file is part of RcppGSL.
##
## RcppGSL is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppGSL is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

NAMESPACE <- environment()

CxxFlags <- function( print = TRUE ){
	inc <- system.file('include', package="RcppRandomSFMT")
	flags <- paste("-I", inc, sep='')
    if( print) cat( flags ) else flags
}

#LdFlags <- function( print = TRUE ){
#	lib <- system.file('libs', package="RcppRandomSFMT")
#	lflags <- paste("-L", lib, sep='')
#	rpath <- paste("-Wl,-rpath,", lib, sep='')
#	flags <- paste(lflags, "-lRcppRandomSFMT", rpath)
#	if( print) cat( flags ) else flags
#}