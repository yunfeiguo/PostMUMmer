#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw/basename/;

die "Usage: $0 <bedtools> <bed> <genome length file> <genome fasta>\n" unless @ARGV == 4;
my $bedtools = shift;
my $bed = shift;
my $genome = shift;
my $fa = shift;
my $outFa = basename $bed.".complement.fasta";
my $complementBED = basename $bed.".complement.bed";

!system("$bedtools sort -i $bed | $bedtools complement -i $bed -g $genome > $complementBED") or die "bedtools complement: $!\n";
!system("$bedtools getfasta -fi $fa -bed $complementBED -fo $outFa") or die "bedtools getfasta: $!\n";
