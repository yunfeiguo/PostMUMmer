#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw/basename/;
use File::Spec;

my $chr_dir = "/home/yunfeiguo/database/hg_index/ucsc/hg38_chr/chroms";
my $tmpdir = "/tmp";
my $pwd = $ENV{PWD};
my $query = File::Spec->catfile($pwd,"hx1.20150513.fasta");
my $qsub_option = <<SCRIPT;
#!/bin/bash
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
#\$ -l h_vmem=8G
#\$ -m a
#\$ -M guoyunfei1989\@gmail.com
set -e
SCRIPT
#\$ -m ea
#\$ -M guoyunfei1989\@gmail.com

for my $i(glob File::Spec->catfile($chr_dir,"*.fa")) {
    #my ($prefix) = (basename $i)=~/^(chr[\dXYM]{1,2})\.fa$/ or warn "ERROR: failed to match $i\n" and next;
    my ($prefix) = (basename $i)=~/^(chr.*?)\.fa$/ or warn "ERROR: failed to match $i\n" and next;
    $prefix = "hx1_$prefix";
    my $script = File::Spec->catfile($tmpdir,"$prefix.".rand($$).".runmummer.sh");
    !system("mkdir $prefix") or die "mkdir(): $!\n" unless -d $prefix;
    open OUT,">",$script or die "ERROR: failed to write to $script: $!\n";
    print OUT $qsub_option,"\n";
    #redirect stderr and stdout
    print OUT "#\$ -e ".File::Spec->catfile($prefix,"stderr")."\n";
    print OUT "#\$ -o ".File::Spec->catfile($prefix,"stdout")."\n";
    #use file staging to reduce IO
    print OUT "cd \$TMP\n";
    print OUT "nucmer --prefix=$prefix $i $query\n";
    print OUT "show-coords -r -c -l $prefix.delta > $prefix.coords\n";
    print OUT "show-tiling $prefix.delta > $prefix.tiling\n";
    print OUT "if [ -s $prefix.tiling ]; then mummerplot --postscript $prefix.tiling -p $prefix;fi\n";
    #print OUT "delta-filter -m $prefix.delta > $prefix.filtered.delta\n";
    #print OUT "show-coords -r -c -l $prefix.filtered.delta > $prefix.filtered.coords\n";
    #print OUT "show-tiling $prefix.filtered.delta > $prefix.filtered.tiling\n";
    #print OUT "if [ -s $prefix.filtered.tiling ]; then mummerplot --postscript $prefix.filtered.tiling -p $prefix.filtered;fi\n";
    print OUT "cp -r * ".File::Spec->catfile($pwd,$prefix)." \n";
    print OUT "touch ".File::Spec->catfile($pwd,$prefix,"$prefix.done")."\n";
    close OUT;
    !system("qsub $script\n") or die "$script: $!\n" unless -e File::Spec->catfile($pwd,$prefix,"$prefix.done");
    #warn("qsub $script\n") unless -e File::Spec->catfile($pwd,$prefix,"$prefix.done");
}
warn "All done\n";


#process coords file
#perl -ne 'next if /^[=\/\s]|^NUCMER/;s/\|//g;print' hx1_50kb_on_hg38.filtered.coords > hx1_50kb_on_hg38.filtered.whitespace.coords
