#! /usr/bin/env perl

# generate the SFMT classes for all combinations of Exponent and for all implementations

use strict;
use warnings;

my @exponents = (607, 1279, 2281, 4253, 11213, 19937, 44497, 86243, 132049, 216091);
my %implementations = (
    NORMAL => { when => 'ALWAYS', define => 'HAVE_NORMAL' },
    SSE2 => { when => 'HAVE_SSE2',   define => 'HAVE_SSE2' },
    ALTIVEC => { when => 'HAVE_AVX', define => 'HAVE_ALTIVEC' }
);


my @headers;
my @sources;
foreach my $mexp (@exponents) {
    foreach my $impl (keys %implementations) {
        #my $class_name = "Sfmt_${mexp}_$impl";
        my $namespace = "SFMT_${impl}_$mexp";

        my $header = "$namespace.hpp";
        my $source = "$namespace.cpp";
        push @headers, $header;
        push @sources, $source;

        my $when = $implementations{$impl}->{when};
        my $define = $implementations{$impl}->{define};

        open(my $fh, ">", $header ) or die "unable to create file [$header]: $!";
        warn "writing file $header\n";

        print $fh <<EOF;
/* This file was automatically generated using perl script $0. Do not edit ! */
#ifndef ${namespace}_HPP
#define ${namespace}_HPP

#include <inttypes.h>
namespace $namespace {
#include "RNG.hpp.t"
}
#endif
EOF
        close $fh;

        open($fh, ">", $source ) or die "unable to create file [$source]: $!";
        warn "writing file $source\n";

        print $fh <<EOF;
/* This file was automatically generated using perl script $0. Do not edit ! */

#include "config.h"
#define ALWAYS

#ifdef $when

#undef HAVE_SSE2
#undef HAVE_ALTIVEC

#define $define
#define MEXP $mexp

//#include "$header"
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <string>
#include <inttypes.h>
namespace $namespace {
#include "RNG.cpp.t"
}
#endif
EOF
        close $fh;

    }
}

my $all_sfmt = "SFMT_include_all.hpp";
open(my $fh, ">", $all_sfmt ) or die "unable to create file [$all_sfmt]: $!";
print $fh <<EOF;
/* This file was automatically generated using perl script $0. Do not edit ! */
#ifndef SFMT_INCLUDE_ALL_HPP
#define SFMT_INCLUDE_ALL_HPP

EOF

print $fh qq{#include "$_"\n} foreach @headers;
print $fh qq{#endif\n};


close $fh;
