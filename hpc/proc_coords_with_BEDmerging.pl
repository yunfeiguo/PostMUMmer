#!/usr/bin/env perl

use strict;
use warnings;
use lib "/home/rcf-02/yunfeigu/perl5";
use SeqMule::Utils;
use SeqMule::Parallel;
use Getopt::Std;
use Data::Dumper;

die "Usage: $0 [options] <split|proc> <indexed query FASTA> [1.coords 2.coords ...]\n".
" -c <TEXT>	contig name\n" 
unless @ARGV >= 2;
#PARAMETERS
my $gap = 100; #max gap allowed for two alignments to be stiched together
my $tmpdir = "/tmp";
my $bedtools = "/home/rcf-02/yunfeigu/proj_dir/Downloads/bedtools2/bin/bedtools";
my $debug = 0;
my $minIdt = 90; #mininum identity between two sequences in a mapping, in percentage
my $minLen2 = 50; #min length of alignment
my $maxOverlapDist = 100; #max distance between a mapping and a CNV for them to be considered as overlapping
my $qsub = "qsub -S /bin/bash -V -l walltime=1:0:0 -l nodes=1:ppn=1 -l mem=4GB -A lc_kw -q laird";
my $cwd = $ENV{PWD};

#INPUT
my %options;
getopts("c:",\%options);
print "\@ARGV: @ARGV\n" if $debug;
my $operation = shift @ARGV;
my $query = shift @ARGV;
my $idx = "$query.fai";
die "no index: $idx\n" unless -e $idx;
my %fa = &SeqMule::Utils::readFastaIdx($idx);
my $total = scalar keys %fa;
my @cleanQ; #files to be removed

if ($operation eq 'split') {
	#split the input coords by contig ID, then submit processing request by qsub
	my $dir = "coord_by_contig";
	mkdir $dir unless -d $dir;
	&splitCoord({fa=>\%fa,dir=>$dir,coord=>\@ARGV});
	for my $i(keys %fa) {
		my $coord = $fa{$i}->{coord} || "";
		my $stderr = File::Spec->catfile($dir,"$i.stderr");
		my $stdout = File::Spec->catfile($dir,"$i.stdout");
		my $mr = File::Spec->catfile($dir,"$i.mr"); #mapping ratio result
		my $done = File::Spec->catfile($dir,"$i.done");
		next if -f $done;
		!system("$qsub -e $stderr -o $stdout ".&SeqMule::Parallel::genTempScript(
				"cd $cwd",
				"$0 -c $i proc $query $coord",
				"touch $done")) or die "qsub $i: $!\n";
		print "we got a coord for this contig\n" if $coord and $debug;
		print("$qsub -e $stderr -o $stdout ".&SeqMule::Parallel::genTempScript(
				"cd $cwd",
				"$0 -c $i proc $query $coord > $mr",
				"touch $done"),"\n") or die "qsub $i: $!\n" if $debug;
	}
} elsif ($operation eq 'proc') {
	my $outputDir = "mapped_query_bed";
	mkdir $outputDir unless -d $outputDir;
	print "\@ARGV: @ARGV\n" if $debug;
	my $coord = shift @ARGV;
	&readCoord(\%fa,$coord);
	print Dumper(%fa) if $debug;
	warn "reading coords done\n";
	&convert2BED(\%fa,$options{c});
	warn "conversion to bed done\n";
	&mergeBED_extractMaxMapping(\%fa,$options{c});
	warn "merging and extraction done\n";
	&output(\%fa,$options{c},$outputDir);
} else {
	die "$operation unknown, use split or proc\n";
}
warn "Clean up ...\n";
&cleanup();
	warn "All done\n";

###################################################################
sub splitCoord {
	#read all coords files, split them into coords with single contig
	my $opt = shift;
	my $fa = $opt->{fa};
	my $coord = $opt->{coord};
	my $dir = $opt->{dir};
	my %fh; #store open filehandles
	my @return; #return individual coords files

	for my $i(@$coord) {
		my %parsedContigs;
		&readCoord(\%parsedContigs,$i);
		for my $j (keys %parsedContigs) {
			unless (defined $fh{$j}) {
				my $coord = File::Spec->catfile($dir,$j.".coords");
				$fa->{$j}->{coord} = $coord;
				#because we are going to append, we should make sure we begin with an empty file
				unlink $coord or die "unlink($coord): $!\n" if -f $coord;
				open ($fh{$j}, ">>",$coord) or die "$coord: $!\n";
			}
			if(defined $parsedContigs{$j}->{raw}) {
				print {$fh{$j}} $_,"\n" for @{$parsedContigs{$j}->{raw}};
			}
		}
	}
	close $_ for values(%fh);
}
sub readCoord {
	#read coords file, store parsed result in a hash
	my $fa = shift;
	my $coord = shift;
	print "will parse $coord\n" if $debug;
	return unless $coord;
	open IN,'<',$coord or die "open($coord): $!\n";
	while(<IN>){
		s/\|//g;
		s/^[\t ]+//;
		next if /^[=\/\s\[]|^NUCMER/;
		chomp;
		my $parsedLine = &parseLine($_);
		my $id = $parsedLine->{id};
		next unless &passFilter($parsedLine);
		($parsedLine->{query_start},$parsedLine->{query_end}) = &smallerFirst($parsedLine->{query_start},$parsedLine->{query_end});
		($parsedLine->{ref_start},$parsedLine->{ref_end}) = &smallerFirst($parsedLine->{ref_start},$parsedLine->{ref_end});
		$parsedLine->{query_start} -= 1; #convert to 0-based start
		$parsedLine->{ref_start} -= 1;
		#we need to make sure all regions that are accepted as final mapped regions come from
		#a continuous sequence on the chr
		if(defined $fa->{$id}->{queryMapping} ) {
			push @{$fa->{$id}->{refMapping}},[$parsedLine->{chr},$parsedLine->{ref_start},$parsedLine->{ref_end}];
			push @{$fa->{$id}->{queryMapping}},[$parsedLine->{chr},$parsedLine->{query_start},$parsedLine->{query_end}];
		} else {
			$fa->{$id}->{refMapping} = [[$parsedLine->{chr},$parsedLine->{ref_start},$parsedLine->{ref_end}]];
			$fa->{$id}->{queryMapping} = [[$parsedLine->{chr},$parsedLine->{query_start},$parsedLine->{query_end}]];
		}
		#store original data, only for splitting
		if (defined $fa->{$id}->{raw}) {
			push @{$fa->{$id}->{raw}},$_;
		} else {
			$fa->{$id}->{raw} = [$_];
		}
		print "Parsing line $.: $_\n" if $debug;
		print Dumper($parsedLine) if $debug;
	}
	close IN;
}
sub passFilter {
	#filter line by global filters
	my $parsedLine = shift;
	my $result = 1;
	$result = 0 if $parsedLine->{idt} < $minIdt;
	$result = 0 if $parsedLine->{len2} < $minLen2;

	return $result;
}
sub parseLine {
	#parse a line in *.coords
	my $line = shift;
	my $result = {};
	my @f=split ' ',$line;
	die "ERROR: expected 13 fields: $line\n" unless @f==13;
	warn "DEBUG:@f\n" if $debug;
#    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
#===============================================================================================================================
#15007271 15007689      32423    32005        419      419      86.64   133797422    51437       0.00     0.81   chr10	m150320_100742_42199_c100794652550000001823158109091525_s1_p0/54078/7931_59368
#0		1	2	3		4	5	6	7		8	  9	  10	  11	12
	#here we only care how much of the query is aligned, regardless of alignment location
	$result->{query_start} = $f[2];
	$result->{query_end} = $f[3];
	$result->{ref_start} = $f[0];
	$result->{ref_end} = $f[1];
	$result->{id} = $f[12];
	$result->{chr} = $f[11];
	$result->{idt} = $f[6];
	$result->{len2} = $f[5]; #alignment length on query

	return $result;
}
sub convert2BED {
	#convert alignment into BED format
	my $ref = shift;
	my $contig = shift;
	my $count = 0;
	for my $i(keys %$ref) {
		next unless $i eq $contig;
		$count++;
		warn "conversion:$count/$total done\n";
		next unless defined $ref->{$i}->{refMapping};
		#the mapping may contain overlapping regions
		#we should only use non-overlapped regions
		my $ref_bed = "$tmpdir/$i".rand($$).".ref.bed";
		my $ref_overlapBED = &SeqMule::Utils::genBED($ref->{$i}->{refMapping});
		$ref->{$i}->{refBED} = &mergeBED($ref_overlapBED,0,$ref_bed);
		push @cleanQ,$ref_overlapBED,$ref_bed;
	}
}
sub mergeBED_extractMaxMapping {
	my $ref = shift;
	my $contig = shift;
	my $count = 0;
	for my $i(keys %$ref) {
		next unless $i eq $contig;
		$count++;
		warn "merge and extraction: $count/$total done\n";
		next unless defined $ref->{$i}->{refBED};
		&findMapping($ref->{$i});
	}
}
sub findMapping {
	my $query_ref = shift;
	my $ref_bed = $query_ref->{refBED};
	my $out = &mergeBED($ref_bed,$gap); #merge with gap allowance
	push @cleanQ,$out;
	#maxMapping = [chr,start,end]
	$query_ref->{maxMapping} = &findMax([&SeqMule::Utils::readBED($out)]);
	($query_ref->{mappedLen},$query_ref->{mappedQueryBed}) = &getMappedLen($query_ref->{maxMapping},$query_ref->{refMapping},$query_ref->{queryMapping});
}
sub getMappedLen {
	#based on maxMapping found on ref
	#identify mapped regions on query
	#get total length of mapped regions on
	#query
	my $max = shift;
	my $allRef = shift;
	my $allQuery = shift;
	my $chr = $max->[0];
	my $start = $max->[1];
	my $end = $max->[2];
	my $mappedLen = 0;
	my @mappedQuery;
	for my $i(0..$#{$allRef}) {
		#the mapped region must be fully contained
		#as we have merged all mapped regions
		next unless $allRef->[$i]->[0] eq $chr and $allRef->[$i]->[1] >= $start and $allRef->[$i]->[2] <= $end;
		push @mappedQuery,$allQuery->[$i];
	}
	#mapped regions in query should be merged with gap of 0 before
	#length is calculated
	my $query_mapped_bed = &SeqMule::Utils::genBED(\@mappedQuery);
	my $query_mapped_nonoverlap_bed = &mergeBED($query_mapped_bed);
	$mappedLen = &SeqMule::Utils::bed2total($query_mapped_nonoverlap_bed);

	push @cleanQ,$query_mapped_bed,$query_mapped_nonoverlap_bed;
	return ($mappedLen,$query_mapped_nonoverlap_bed);
}
sub findMax {
	#given a set of regions in BED format
	#find the longest the region
	my $ref = shift;
	my $max = 0;
	my $pos = [];
	#array of [chr,start,end]
	for my $i(@$ref) {
		my $len = $i->[2]-$i->[1];
		#if there are multiple alignments with equal lengths,
		#we output first one of them
		if ($len > $max) {
			$pos = $i;
			$max = $len;
		}
	}
	return $pos;
}
sub output {
	my $ref = shift;
	my $contig = shift;
	my $outputDir = shift;
	warn "mapped regions in query will be saved to $outputDir in BED format.\n";
	my %fa = %$ref;
	print join("\t","FASTA_ID","Mapped_length","Total_length","Mapping_ratio","Mapped_chr","Query_start","Query_end"),"\n";
	for my $i(keys %fa) {
		next unless $i eq $contig;
		my $mappedLen = $fa{$i}->{mappedLen} || 0;
		my $totalLen = $fa{$i}->{length};
		my $mapping = ['NA','NA','NA'];
		if(defined $fa{$i}->{maxMapping}) {
			$mapping = $fa{$i}->{maxMapping};
		}
		my $mappingRatio = $totalLen == 0? 0:$mappedLen/$totalLen;
		#output mapping stats
		print join("\t",$i,
			$mappedLen,
			$totalLen,
			sprintf("%.2f",$mappingRatio),
			@$mapping,
		), "\n";
		#write mapped regions of query in BED format
		#we need to replace chr name by query contig name
		#because previously we put chr name here to differentiate 
		#alignments on different chromosomes.
		!system("perl -pe 's/^(\\S+)/$i/' ".$fa{$i}->{mappedQueryBed}." > $outputDir/$i.mappedQuery.bed")
		 or die "failed to copy BED $i: $!\n" if defined $fa{$i}->{mappedQueryBed};
	}
}
sub cleanup {
	#remove temp files
	unlink @cleanQ;
}
sub smallerFirst {
	my $start = shift;
	my $end = shift;
	if($start > $end) {
		#make sure end is no smaller than start
		$start = $end + $start;
		$end = $start - $end;
		$start = $start - $end;
	}
	return($start,$end);
}
sub mergeBED {
	my $bed = shift;
	my $gap = shift;
	my $out = shift;
	$gap = $gap || 0;
	$out = $out || "$tmpdir/$$".rand($$).".tmp.bed";
	!system("$bedtools sort -i $bed | $bedtools merge -d $gap > $out") or die "merging $bed fail: $!\n";
	return $out;
}
