#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# Usage
my $usage = "
fq2fa.pl : converts fastq to fasta
Usage: perl fq2fa.pl -i <fastq> -o <fasta>
";
my($inf, $outf);
GetOptions(
    'i:s'         => \$inf,
    'o:s'         => \$outf,
);
die $usage unless ($inf);

my ($temp, $in, $out);
my @temp;

# Open input and output files.
open( $in, "<",  $inf)  or die "Can't open $inf: $!";
open( $out, ">",  $outf)  or die "Can't open $outf: $!";

while (<$in>){
  chomp($temp[0] = $_);		# First line is an id.
  chomp($temp[1] = <$in>);	# Second line is a sequence.
  chomp($temp[2] = <$in>);	# Third line is an id.
  chomp($temp[3] = <$in>);	# Fourth line is quality.

  # Prune first char.
  $temp[0] = substr($temp[0], 1);

  # Substring to inset value.
 # $temp[1] = substr($temp[1], $inset);

  # Print to fasta file.
  print $out ">$temp[0]\n";
  print $out "$temp[1]\n";
}

close $in or die "$in: $!";
close $out or die "$out: $!";
