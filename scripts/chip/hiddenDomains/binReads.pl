#!/usr/bin/perl

##========================================================================
##
## Author: Joshua Starmer <josh.starmer@gmail.com>, 2014
## 
## Copyright (C) 2014, Joshua Starmer
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
## 
##========================================================================


use strict;
use warnings;

use Getopt::Std;

use constant TRUE  => 1;
use constant FALSE => 0;

use constant SAM_CHR  => 2;
use constant SAM_POS  => 3;
use constant SAM_QUAL => 4;
use constant SAM_SEQ  => 9;

use constant BED_CHR => 0;
use constant BED_POS => 1;

my $sHelp = <<END_OF_HELP;

Usage: binReads.pl [options] [-B bedFile.bed | bamFile.bam] > binned_reads.txt

Options

-h
    Print this help information.

-b BIN_WIDTH
    The width of the bin. Default is 1000bp.

-B
    The input file is in BED format (the default is BAM)

-q  MIN_MAPQ
    The minimum MAPQ score. Default is 30.

-M
    Assume all bins should be on mouse chromosomes. This is the default.

-H
    Assume all bins should be on human chromosomes.

-c  "chr1 chr2 ..."
    Bin reads only from specified chromosomes.
    

END_OF_HELP

if (-t STDIN && !@ARGV) {
    print STDERR $sHelp;
    exit 1;
}


# process the command line arguments
my %opts; # a hash table to store file names passed in as agruments
getopts('hb:mMHc:q:B', \%opts);

if ($opts{'h'}) { # print help and exit
    print STDERR $sHelp;
    exit 1;
}

my $binWidth = 1000;
if (defined($opts{'b'})) {
    $binWidth = $opts{'b'};
}
print STDERR "binWidth: ".$binWidth." (change with -b option)\n";

my $minQualScore = 30;
if (defined($opts{'q'})) {
    $minQualScore = $opts{'q'};
}
print STDERR "minQualScore: ".$minQualScore." (change with -q option)\n";


my @chrList = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY");

if (defined($opts{'m'}) || defined($opts{'M'})) {
    # use mouse chromosomes
    print STDERR "Using mouse chromosomes\n";
    @chrList = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY");
} elsif (defined($opts{"H"})) {
    print STDERR "Using human chromosomes\n";
    @chrList = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY");
} elsif (defined($opts{"c"})) {
    print STDERR "Using custom chromosome list\n";
    @chrList = split(" ", $opts{"c"});
} else {
    print STDERR "Default: Using mouse chromosomes.  Change this with -m, -h or -c\n";
}


my $bamFile = shift(@ARGV);
my $INPUT_FILE;

my $usingBedFile = FALSE;
if(defined($opts{'B'})) {
    $usingBedFile = TRUE;
    open($INPUT_FILE, "<".$bamFile) || die("Could not open $bamFile $!\n");
} else {
    open($INPUT_FILE, "samtools view -F 4 -q $minQualScore $bamFile |") 
	|| die("Could not open $bamFile $!\n");
}
my %chrBins;
my $counter = 0;
print STDERR "Binning reads.";
while(my $alignment = <$INPUT_FILE>) {
    if ($alignment =~ /^(\@)/) { # skipp header lines...
	next;
    }

    if (($counter % 100000) == 0) {
	print STDERR ".";
    }
    
    chomp($alignment);

    my $chr;
    my $pos;
    if ($usingBedFile) {
	my @values = split(/\t/, $alignment);
	
	$chr = $values[BED_CHR];
	$pos = $values[BED_POS];
    } else {
	my @sam = split(/\t+/, $alignment); 

	$chr = $sam[SAM_CHR];
	$pos = $sam[SAM_POS] - 1;
#    my $quality = $sam[SAM_QUAL];
#    my $seqLen = length($sam[SAM_SEQ]);
    }

    my $binIndex = sprintf("%d", $pos / $binWidth) * $binWidth;    

    if (defined($chrBins{$chr})) {
	my $binsRef = $chrBins{$chr};
	if (defined($$binsRef{$binIndex})) {
	    $$binsRef{$binIndex}++;
	} else {
	    $$binsRef{$binIndex} = 1;
	}
    } else {
	my %bins;
	$bins{$binIndex} = 1;
	$chrBins{$chr} = \%bins;
    }

    $counter++;
}


print STDERR "\nbinned ".$counter." reads\n";
if ($counter == 0 && !$usingBedFile) {
    print STDERR "\nNone of the reads in the BAM file passed the quality";
    print STDERR " threshold, ".$minQualScore.".\nUse the -q option to set ";
    print STDERR "a lower threshold\n\n";
    exit;
}

print "id\tchr\tpos\tcount\n";
$counter = 0;
foreach my $chr (@chrList) {
    #print "chr: ".$chr."\n";
    my $binsRef = $chrBins{$chr};
    foreach my $binIndex (sort {$a <=> $b} keys(%$binsRef)) {
	print $counter."\t".$chr."\t".$binIndex."\t".$$binsRef{$binIndex}."\n";
	$counter++;
    }
}

