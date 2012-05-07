This package wraps the original SFMT implementation by Mutsuo Saito and Makoto Matsumoto
available from <http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/>.

It makes some changes because the original implementation uses global state for the RNG,
that makes it impossible to use for instance in multithreaded code.

Here are some instructions on how the original implementation was wrapped, that may be useful
to wrap an updated version:

1) all the files from the original implementation are in src/SFMT-src

2) the only source file is SFMT.hpp
Look in this files for comments "WRAP INSTRUCTIONS" that explain where and how to copy code

The originals file form which you will have to copy/paste code are:
    - SFMT.c
    - SFMT-sse2.h
    - SFMT-alti.h
