#! /usr/bin/env perl

# generate the SFMT classes for all combinations of Exponent and for all implementations

use strict;
use warnings;

my $AUTO = 'auto';
chdir $AUTO or die "impossible to set current directory to $AUTO";

my @exponents = (607, 1279, 2281, 4253, 11213, 19937, 44497, 86243, 132049, 216091);
my %implementations = (
    NORMAL => { when => 'ALWAYS', define => 'HAVE_NORMAL' },
    SSE2 => { when => 'HAVE_SSE2',   define => 'HAVE_SSE2' },
    ALTIVEC => { when => 'HAVE_AVX', define => 'HAVE_ALTIVEC' }
);

my @undefs = qw(MEXP SFMT_H SFMT_PARAMS_H SFMT_ALTI_H SFMT_SSE2_H N N32 N64
    POS1 SL1 SL2 SR1 SR2 MSK1 MSK2 MSK3 MSK4 PARITY1 PARITY2 PARITY3 PARITY4
    ALTI_SL1 ALTI_SR1 ALTI_MSK ALTI_MSK64 ALTI_SL2_PERM ALTI_SL2_PERM64 ALTI_SR2_PERM ALTI_SR2_PERM64
    IDSTR
);
#undef MEXP
#undef SFMT_CLASS_NAME
#undef SFMT_H
#undef SFMT_PARAMS_H

#undef SFMT_ALTI_H
#undef SFMT_SSE2_H
my @files;
foreach my $mexp (@exponents) {
    foreach my $impl (keys %implementations) {
        #my $class_name = "Sfmt_${mexp}_$impl";
        my $namespace = "SFMT_${mexp}_$impl";

        my $class_name = "RNG";
        my $file = "Sfmt_${mexp}_$impl.hpp";
        push @files, $file;
        open(my $fh, ">", $file ) or die "unable to create file [$file]: $!";
        warn "writing file $file\n";

        my $when = $implementations{$impl}->{when};
        my $define = $implementations{$impl}->{define};

        print $fh <<EOF;
/* This file was automatically generated using perl script $0. Do not edit ! */
#ifndef ${namespace}_HPP
#define ${namespace}_HPP

#define ALWAYS
#ifdef $when // otherwise not supported

// force neutral implementation
#ifdef HAVE_SSE2
#define RESET_HAVE_SSE2
#undef HAVE_SSE2
#endif

#ifdef HAVE_ALTIVEC
#define RESET_HAVE_ALTIVEC
#undef HAVE_ALTIVEC
#endif

#define $define 1

#ifndef READ_CONFIG_H
#include "../config.h"
#define READ_CONFIG_H
#endif

#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#ifdef HAVE_EMMINTRIN_H
#include <emmintrin.h>
#endif

#ifdef HAVE_ALTIVEC_H
#include <altivec.h>
#endif

#define MEXP $mexp


#include <string>
#include <string.h>
#include <assert.h>
namespace $namespace {
#include "../Sfmt_base.hpp"
}

#undef $define

#ifdef RESET_HAVE_SSE2
#define HAVE_SSE2
#endif

#ifdef RESET_HAVE_ALTIVEC
#define HAVE_ALTIVEC_H
#endif

#undef SFMT_PARAMS${mexp}_H

EOF

        print $fh "#undef $_\n" foreach (@undefs);
        print $fh "#endif // when\n";
        print $fh "#endif // PROTECT\n"; # last endif


        close $fh;
    }
}

my $all_sfmt = "Sfmt_include_all.hpp";
open(my $fh, ">", $all_sfmt ) or die "unable to create file [$all_sfmt]: $!";
print $fh <<EOF;
/* This file was automatically generated using perl script $0. Do not edit ! */
#ifndef SFMT_INCLUDE_ALL_HPP
#define SFMT_INCLUDE_ALL_HPP
#include <inttypes.h>
EOF

print $fh qq{#include "$_"\n} foreach @files;
print $fh qq{#endif\n};


close $fh;
