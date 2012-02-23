#! /usr/bin/env Rscript
#############################################################################
# run_runit.R
# Run the Runit tests of a package either via R CMD check or by directly running this script
# It is greatly inspired by doRUnit.R from http://rwiki.sciviews.org/doku.php?id=developers:runit
#
### Running the tests:
# 1) via R CMD check
#	 this script has to be stored in the tests/ directory. It will find the runit tests and execute them
#    and produce a fatal error if any test fail (stop()) to interrupt the check process
#
# 2) running the script directly in a source package directory (uninstalled)
#    Because the check process is very slow, especially for test-driven development
#    we developed another way to run the tests.
#    In this process, your source files will be sourced (in the Collate order if specified)
#    then the unit tests will be executed.
#    If an .onLoad() function is defined in the sources, it will be executed prior
#    to run the tests.
#    WARNING: this is not equivalent to the R CMD check because:
#      - the package is never installed, even temprorarily
#      - so it can not work with foreign code
#      - so the namespace, exports, imports etc.. are not done
#
# 3) running the script directly in a installed package directory (in library/)
#    You can re-run the test suite from an installed package. It will load the library
#    then run the unit tests.



### Setup

# 1) copy this script in the tests/ directory of your package
# 2) put your runit scripts in inst/unitTests/ with filenames beginning with runit.
# 3) copy this script in inst/tests/ to be sure it is installed by R CMD INSTALL even if the --install-tests
#    option is not used



DEBUG <- TRUE

# should be true iff sourced from a console
SOURCED <- length(sys.frames()) > 0
main <- function() {
	## unit tests will not be done if RUnit is not available
	if( ! require("RUnit", quietly=TRUE)) {
		warning("cannot run unit tests -- package RUnit is not available")
		return()
	}

	source_dir <- ""
	installed_loc <- ""
	runit_dir <- ""

	rscript_path <- findRscriptPath()
	rscript <- nchar(rscript_path) > 0

	if ( rscript ) {
		# we are in a direct use with optionally  an argument
		# it could in a source package or in a already installed package (e.g. post-installation tests)

		pkg_dir <- parentDir( dirname(rscript_path) )
		if ( pkg_dir == "." || pkg_dir == "")
			pkg_dir = getwd()

		args <- commandArgs(TRUE)
		if ( length(args) ) {
			installed_loc <- args[1]
			source_dir <- pkg_dir
		} else {
			if ( isInstalledPackageDirectory(pkg_dir) ) {
				installed_loc <-parentDir(pkg_dir)
				source_dir <- ""
			} else {
				source_dir <- pkg_dir
				installed_loc <- ""
			}
		}

	} else { # no Rscript, so we are sourced
		# could be via R CMD check or by a manual source() in a console
		pkg_dir <- findPackageDirectory()

		if ( SOURCED ) {
			source_dir <- pkg_dir
			installed_loc <- ""
		} else { # R CMD check
			cat("...R CMD CHECK CONTEXT DETECTED\n")
			source_dir <- ""
			installed_loc <- parentDir(pkg_dir)
			runit_dir <- file.path(pkg_dir, "unitTests")
			# SPECIAL CASE for R CMD check and data files location
			cat("...SETTING CURRENT DIRECTORY TO ", pkg_dir, "\n")
			setwd(pkg_dir)


		}

	}

	paths <- strsplit(pkg_dir, "/")[[1]]
	pkg <- paths[length(paths)]

	if (SOURCED){
		RUNIT_TRUE_TEST <<- function() {
			runit(pkg, source_dir = source_dir, installed_loc = "..truetest")
		}
	}

	if (DEBUG) {
		cat('----------- DEBUG -------------------\n')
		cat("getwd()=", getwd(), "\n")
		cat("Guessed package directory = ", pkg_dir, "\n")
		cat("rscript_path =", rscript_path, "\n")
		cat("Are we in a Rscript context ? =", rscript, "\n")
		cat("source_dir=", source_dir, "\n")
		cat("installed_loc=", installed_loc, "\n")
		cat("runit_dir=", runit_dir, "\n")
		cat("pkg=", pkg, "\n")
	}

	runit(pkg, source_dir = source_dir, installed_loc = installed_loc, runit_dir = runit_dir)


}

# pkgname: the name of the package
# source_dir: where the sources are (sometimes not available)
runit <- function(pkgname, source_dir = "", installed_loc = "", runit_dir = "") {
	if ( ! nchar(source_dir) && ! nchar(installed_loc) )
		stop("Bad params, need at least one source or one installed_loc")

	if ( pkgname %in% loadedNamespaces() ) {
		cat("unloading the namespace\n")
		unloadNamespace(pkgname)
	}


	### loading the package ###
	if ( nchar(installed_loc)  ) {
		# we have an installed library to use, so let's load it
		cat("loading library ", pkgname, " from lib.loc=", installed_loc, "\n")
		library(package=pkgname, character.only=TRUE, lib.loc = installed_loc)
	} else {
		if ( DEBUG ) {
			cat("====================================================================\n")
			cat("loading SOURCE package in ", source_dir, "\n")
			cat("====================================================================\n")
		}
		loadSourcePackage(pkgname, source_dir)
	}


	### test suite
	# guess runit_dir
	if ( ! nchar(runit_dir) ) {
		if ( nchar(source_dir) ) {
			runit_dir <- file.path(source_dir, "inst", "unitTests")
		} else {
			runit_dir <- file.path(installed_loc, "unitTests")
		}
	}

	# Define tests
	testSuite <- defineTestSuite(name=paste(pkgname, "runit testing"), dirs=runit_dir)
	# Run it
	tests <- runTestSuite(testSuite)
	# Report
	cat("------------------- UNIT TEST RESULTS ---------------------\n\n")
	printTextProtocol(tests, showDetails=TRUE)
	cat("\n\n------------------- UNIT TEST SUMMARY ---------------------\n\n")
	printTextProtocol(tests, showDetails=FALSE)

	## Return stop() to cause R CMD check stop in case of
	##  - failures i.e. FALSE to unit tests or
	##  - errors i.e. R errors
	tmp <- getErrors(tests)
	if(tmp$nFail > 0 | tmp$nErr > 0) {
		stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
						", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
	}


}

## Return true iff the current script is executed via Rscript
## current heuristic is to detected the --args option used by Rscript
#inRscript <- function() {
#	return( length( grep("--args", commandArgs() ) ) > 0 )
#}

# get parent dir
parentDir <- function(path) {
	paths <- strsplit(path, "/")[[1]]
	nb <- length(paths)
	return( paste(paths[1:nb-1], collapse="/") )
}

# return the path of the Rscript script if executed via Rscript
# see http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
findRscriptPath <- function() {
	initial.options <- commandArgs(trailingOnly = FALSE)
	file.arg.name <- "--file="
	script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
	if ( length(script.name) == 0)
		script.name <- ""
	return(script.name)
}

# tells if the pkg directory is an installed one or a source one
isInstalledPackageDirectory <- function(pkg_dir) {
	# heuristic: is there a Meta directory ?
	meta <- file.path(pkg_dir, "Meta")
	return(file.exists(meta))
}


## do all the prepartion needed to run the unit tests
## via the source() method
#prepare_runtests_via_source <- function(pkg_dir) {
#	Rsrc   <- file.path(pkg_dir ,'R')
#	description <- file.path(pkg_dir, "DESCRIPTION")
#
#	####### description file: to get the collate field
#
#	if ( ! file.exists(description) )
#		stop("No DESCRIPTION file found in directory ", basename)
#	if ( DEBUG )
#		cat("...found DESCRIPTION file\n")
#	collate <- read.dcf(file=description, field="Collate")[[1]]
#
#	####### get the list of R source files
#	files_to_source <- NULL
#	if ( ! is.na(collate) ) {
#		if ( DEBUG )
#			cat("...Found collate files in DESCRIPTION file\n")
#
#		f <- gsub("\n", " " , collate) # remove end of lines
#		f <- strsplit(f, " ")[[1]]
#		files_to_source <- gsub("'", "" , f) # remove single quotes
#	} else {
#		if ( DEBUG )
#			cat("No Collate found\n")
#		# get all the files in the src directory
#		files_to_source <- sort(list.files(patt='\\.[Rr]$', Rsrc))
#	}
#	# add the src/ directory to the filename
#	files_to_source <- as.character(sapply(files_to_source, function(f) { return(file.path(Rsrc, f)) }))
#	if (DEBUG)
#		cat("Files to source:", files_to_source, "\n")
#
#	##### source the R files
#
#	# create a namespace, source the files into it and export all
#
#
##	loadings <- sapply(files_to_source, source, echo=FALSE, verbose=FALSE)
#	# if a .onLoad() function has been defined, execute it
#
#
#	if ( exists(".onLoad") && is.function(get(".onLoad")) )
#		.onLoad()
#	if( DEBUG )
#		cat("source files sourced\n")
#
#
#
#}

makeNamespace <-function(name, version = NULL, lib = NULL) {
	impenv <- new.env(parent = .BaseNamespaceEnv, hash = TRUE)
	attr(impenv, "name") <- paste("imports", name, sep = ":")
	env <- new.env(parent = impenv, hash = TRUE)
	name <- as.character(as.name(name))
	version <- as.character(version)
	info <- new.env(hash = TRUE, parent = baseenv())
	assign(".__NAMESPACE__.", info, envir = env)
	assign("spec", c(name = name, version = version),
			envir = info)
	setNamespaceInfo(env, "exports", new.env(hash = TRUE,
					parent = baseenv()))
	setNamespaceInfo(env, "imports", list(base = TRUE))
	setNamespaceInfo(env, "path", file.path(lib, name))
	setNamespaceInfo(env, "dynlibs", NULL)
	setNamespaceInfo(env, "S3methods", matrix(NA_character_,
					0L, 3L))
	assign(".__S3MethodsTable__.", new.env(hash = TRUE,
					parent = baseenv()), envir = env)
	.Internal(registerNamespace(name, env))
	env
}


# load a source package like with library but with limited functionality
loadSourcePackage <- function(pkg, pkg_dir=getwd()) {
	which.lib.loc <- parentDir(pkg_dir)

	nsinfo <- parseNamespaceFile(pkg, which.lib.loc)
#	print(nsinfo)
	dataPath <- file.path(which.lib.loc, pkg, "data")
	env <- makeNamespace(pkg)


	Rsrc   <- file.path(pkg_dir ,'R')
	description <- file.path(pkg_dir, "DESCRIPTION")

	####### description file: to get the collate field

	if ( ! file.exists(description) )
		stop("No DESCRIPTION file found in directory ", basename)
	if ( DEBUG )
		cat("...found DESCRIPTION file\n")
	collate <- read.dcf(file=description, field="Collate")[[1]]

	####### get the list of R source files
	files_to_source <- NULL
	if ( ! is.na(collate) ) {
		if ( DEBUG )
			cat("...Found collate files in DESCRIPTION file\n")

		f <- gsub("\n", " " , collate) # remove end of lines
		f <- strsplit(f, " ")[[1]]
		files_to_source <- gsub("'", "" , f) # remove single quotes
	} else {
		if ( DEBUG )
			cat("No Collate found\n")
		# get all the files in the src directory
		files_to_source <- sort(list.files(patt='\\.[Rr]$', Rsrc))
	}
	# add the src/ directory to the filename
	files_to_source <- as.character(sapply(files_to_source, function(f) { return(file.path(Rsrc, f)) }))
	if (DEBUG)
		cat("Files to source:", files_to_source, "\n")

#	loadings <- sapply(files_to_source, source, echo=FALSE, verbose=FALSE)
	for (file in files_to_source)
		sys.source(file, envir = env)

	# if a .onLoad() function has been defined, execute it

	if ( exists(".onLoad", envir = env) && is.function(get(".onLoad", envir = env)) ) {

		if ( DEBUG )
			cat("Executing .onLoad() method\n")
		get(".onLoad", envir = env)()
	}
	if( DEBUG )
		cat("source files sourced\n")

	## exports
	cat("Exports:", nsinfo$exports, "\n")
	namespaceExport(env, nsinfo$exports)

	## imports
#	imports <- names(getNamespaceInfo(env, "imports"))
	imports <- nsinfo$imports[[1]]

#	cat("Imports:", imports, "\n")
	for (l in imports) {
#		print(l)
		cat("loading library ", l, "\n")
		library(l, character.only = TRUE)
	}

	attachNamespace(env, dataPath = dataPath)

	return(env)
}

# return the absolute package directory path
findPackageDirectory <- function(dir=getwd()) {
	paths <- strsplit(dir, "/")[[1]]
#	cat(paths)
	nb <- length(paths)
	pkgdirs <- NULL
	current <- paths[nb]
	if ( current == "tests") {  # we must probably are in the tests subdirectory
		parent <- paths[nb-1]
		if (length( grep(".Rcheck", parent ) )) {
			# we are must probably called by a R CMD check
			# in a .Rcheck directory
			# try to locate the package directory, which contain DESCRIPTION file
			parent_dir <- paste(paths[1:nb-1], collapse="/")
#			print(parent_dir)
			files <- list.files(parent_dir, recursive=TRUE)
			description_path <- files[grep("^\\D.*/DESCRIPTION", files, perl=TRUE)]
			if ( ! length(description_path) )
				stop("unable to find the DESCRIPTION file somewhere in the .Rcheck directory:", parent_dir)
			pkg_dir_name <- strsplit(description_path, "/")[[1]][1]
			if ( ! length( pkg_dir_name ))
				stop("unable to extract the package directory name from the description file:", description_path)
			pkgdirs <- c(paths[1:(nb-1)], pkg_dir_name)

		} else {
			pkgdirs <- paths[1:(nb-1)]
		}


	} else {
		# heuristic: assume we are in the package directory
		pkgdirs <- paths
	}

	return(paste(pkgdirs, collapse="/"))
}

#findRunitDirectory <- function(pkg_dir) {
#	# simple heuristic: if no inst/ subdir, then
#	# we are in an installed directory, so go directly for unitTests
#	dir <- file.path(pkg_dir, "inst")
#	if ( ! file.exists(dir) )
#		dir <- pkg_dir
#
#	path <- file.path(dir, "unitTests")
#
#	return(path)
#}

main()
