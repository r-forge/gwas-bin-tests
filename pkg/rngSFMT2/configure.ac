# Process this file with autoconf to produce a configure script.
AC_INIT(rngSFMT2, 0.1, karl.forner@gmail.com)

# safety check, probably useless that it is the right src directory
AC_CONFIG_SRCDIR([src/config.h.in])

AC_CONFIG_HEADER(src/config.h)


### ** Precious variables.: setting MEXP
AC_ARG_VAR([MEXP],
           [Mersenne exponent and the period of generated code will be 2MEXP-1. MEXP must be one of 607, 1279, 2281, 4253, 11213, 19937, 44497, 86243, 132049, 216091. ])
# setting default value
if test -z "$MEXP"; then
  MEXP=19937;
fi
AC_DEFINE_UNQUOTED([MEXP], [$MEXP], [The Mersenne exponent])

### forcing not to use SSE2
AC_ARG_ENABLE([sse2],
              [AC_HELP_STRING([--disable-sse2],
                              [Disable sse2])],
              [],
              [enable_sse2=yes])


# Checks for programs.
AC_PROG_CC
AC_PROG_CXX

AC_CHECK_HEADERS([inttypes.h emmintrin.h altivec.h])


# check endianness define BIG_ENDIAN if big endian
AC_C_BIGENDIAN(
    [AC_DEFINE([BIG_ENDIAN], 1,[Define IF IT IS bIG ENDIAN])]
              )

# additional macros

m4_include([m4/m4_ax_check_compile_flag.m4])
m4_include([m4/m4_ax_gcc_x86_cpuid.m4])
m4_include([m4/ax_ext2.m4])

# to detect SIMD extensions:
AX_EXT2

# ALTIVEC SUPPORT
if test "$ax_cv_have_avx_ext" = yes; then
    AC_DEFINE(HAVE_ALTIVEC,1,[Support AVX (Advanced Vector Extensions) instructions])
fi

if test "$ax_cv_have_sse2_ext" = yes; then
    echo "HAVE_SSE2=$HAVE_SSE2"
    HAVE_SSE2=1
    if test "x$enable_sse2" = "xno"; then
        echo "SSE2 supported BUT disabled by user"
        AC_DEFINE(HAVE_SSE2, 0)
    fi
fi


AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
