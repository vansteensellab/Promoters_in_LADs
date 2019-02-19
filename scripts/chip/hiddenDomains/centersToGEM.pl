#!/usr/bin/perl

use strict;
use warnings;

while(my $line = <>) {
    chomp($line);
    my @data = split(/\t/, $line);

    my $chr = $data[0];
    my $pos = $data[1];
    
    $chr =~ s/chr//;
    print $chr.":".$pos."\n";
}
