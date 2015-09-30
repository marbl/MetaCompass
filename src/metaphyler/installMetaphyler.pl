#!/usr/bin/perl

#############################################
#
# Program: Install Metaphyler for MetaCompass.
#
# Author: Bo Liu
#
# Tue Aug  7 16:02:48 EDT 2012
#
#############################################

use strict;
use warnings;
use FindBin qw($Bin);


# format blast database
my $cmd = "formatdb -p F -i $Bin/markers/markers.refseq.dna";
print "$cmd\n";
system($cmd);

system("mkdir $Bin/bin");
my $gcc = "g++ -Wall -W -O2";
my @programs = ("metaphylerClassify", "taxprof");
foreach my $program (@programs) {
    $cmd = "$gcc -o $Bin/bin/$program $Bin/src/$program.cpp";
    print "$cmd\n";
    system($cmd);
}

exit;
