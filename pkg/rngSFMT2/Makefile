
PKG=${shell pwd}
PKGNAME=${shell basename ${PKG}}
SUITE=tests/run_runit.R

# As from R 2.13.0, R_HOME/bin is prepended to the PATH so that references to R or Rscript in the Makefile do make use of the currently running version of R.
R=R
RSCRIPT=Rscript

NAMESPACE = NAMESPACE
DESCRIPTION = DESCRIPTION

TRUETEST_DIR = ..truetest

all: ${NAMESPACE}

# generates the manual in pdf: Rd2.pdf

pdf: Rd2.pdf

pdfclean: 
	rm -rf .Rd2dvi*

Rd2.pdf: ${NAMESPACE} ${DESCRIPTION} ${R_SOURCES}
	rm -f Rd2.pdf
	$R CMD Rd2dvi --pdf --no-preview man/
# roxygen: run roxygen that generates the man/*.Rd documentation files and
# the NAMESPACE file
roxygen: ${NAMESPACE} ${DESCRIPTION}

${NAMESPACE} ${DESCRIPTION}: ${R_SOURCES}
	rm -f man/*.Rd
	echo 'library(roxygen); roxygenize(".",".", overwrite=TRUE, copy.package=FALSE, use.Rd2=TRUE)' \
	| $R --quiet --vanilla
	


# Run unit tests
test: ${NAMESPACE} ${DESCRIPTION}
	${RSCRIPT} ${SUITE}

# run R CMD check
check: ${NAMESPACE} ${DESCRIPTION} clean
	$R CMD check .

## check the documentation, but not the installation and the tests
#checkdoc: ${NAMESPACE} ${DESCRIPTION} clean
#	$R CMD check --no-install --no-tests --check-subdirs=no .

${TRUETEST_DIR}:
	mkdir ${TRUETEST_DIR}

${TRUETEST_DIR}/${PKGNAME}/${SUITE}: ${NAMESPACE} ${DESCRIPTION} ${R_SOURCES} ${TEST_SOURCES}
	MAKEFLAGS="-j 6" $R CMD INSTALL --no-docs --install-tests -l ${TRUETEST_DIR} .

# make a real install inside the ..truetest sub directory, then run the tests
# in this installed package
true-test: ${NAMESPACE} ${DESCRIPTION} ${TRUETEST_DIR} ${TRUETEST_DIR}/${PKGNAME}/${SUITE}
	${RSCRIPT} ${SUITE} ${TRUETEST_DIR}

# fast check: in fact not much faster
fast-check: ${NAMESPACE} ${DESCRIPTION}
	$R CMD check --no-codoc --no-examples --no-vignettes --no-latex .

debug:
	@echo R_SOURCES=${R_SOURCES}
	@echo PKG=${PKG} PKGNAME=${PKGNAME}

make notsoclean:
	rm -rf ${TRUETEST_DIR} ..Rcheck
	rm -f  src/*.o src/*.so  */*~ *~

clean: pdfclean
	rm -rf ${TRUETEST_DIR} ..Rcheck
	rm -f  src/*.o src/*.so src/libsrc/*.o */*~ *~
	rm -f config.log config.status

	
