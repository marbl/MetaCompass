#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);

my $cmd = "";
 
print "\n";
print "# Install MetaCompass\n";
$cmd = "g++ -Wall -W -O2 -o ./bin/extractSeq ./src/utils/extractSeq.cpp";
print "$cmd\n";
system($cmd);

$cmd = "g++ -Wall -W -O2 -o ./bin/formatFASTA ./src/utils/formatFASTA.cpp";
print "$cmd\n";
system($cmd);


$cmd = "g++ -Wall -W -O2 -o ./bin/buildcontig ./src/buildcontig/buildcontig.cpp src/buildcontig/cmdoptions.cpp src/buildcontig/memory.cpp src/buildcontig/procmaps.cpp src/buildcontig/outputfiles.cpp";
print "$cmd\n\n";
system($cmd);

