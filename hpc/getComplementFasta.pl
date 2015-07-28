#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw/basename/;
my $minLen = 10_000; #min length for unmapped region

die "Usage: $0 <bed> <indexed genome fasta>\n" unless @ARGV == 2;
my $bedtools = "bedtools";
my $bed = shift;
my $genome = "/tmp/".rand($$).".genome";
my $fa = shift;
my $fai = "$fa.fai";
my $outFa = basename $bed.".complement.$minLen.fasta";
my $complementBED = basename $bed.".complement.$minLen.bed";

!system("cut -f 1,2 $fai > $genome") or die "cut $fai: $!\n";
!system("$bedtools sort -i $bed | $bedtools complement -i $bed -g $genome | perl -ne '\@f=split;print if \$f[2]-\$f[1] >= $minLen' > $complementBED") or die "bedtools complement: $!\n";
!system("$bedtools getfasta -fi $fa -bed $complementBED -fo $outFa") or die "bedtools getfasta: $!\n";
