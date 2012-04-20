# Process this file with autoconf to produce a configure script.
AC_INIT(rngSFMT2, 0.1, karl.forner@gmail.com)

# safety check, probably useless that it is the right src directory
AC_CONFIG_SRCDIR([src/config.h.in])

AC_CONFIG_HEADER(src/config.h)

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

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
