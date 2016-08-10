#!/usr/bin/perl

#############################################
#
# Program: Pick reference genomes based on
#          metaphyler output
#
# Author: Bo Liu
#
# Fri Aug 10 18:12:01 EDT 2012
#
#############################################

use strict;
use warnings;
use FindBin qw($Bin);


#----------------------------------------#
# read command line options
#----------------------------------------#
my $blast = "";
my $part = 0;
if (scalar @ARGV == 2) {
    ($blast, $part) = @ARGV;
} else {
    Usage();
}
#----------------------------------------#


#----------------------------------------#
# information of reference genomes
#----------------------------------------#
my %seq2tid = ();    # sequence id to taxonomy id
my %tid2seqs = ();   # taxonomy id to sequences
my %tid2name = ();   # tax id to name
my %tid2sp = ();     # tax id to species id
open(FH, "$Bin/../refseq/tid2par.tab") or die("Could not open $Bin/tid2par.tab\n");
foreach my $line (<FH>) {
    chomp $line;
    my ($seq, $tid, $sp, $ge, $name) = split("\t", $line);
    $seq2tid{$seq} = $tid;
    $tid2seqs{$tid}{$seq} = 1;
    $tid2name{$tid} = $name;
    if ($sp eq "NA") { $sp = $tid;}
    $tid2sp{$tid} = $sp;
}
close FH;
#----------------------------------------#


#----------------------------------------#
# abundance for each reference genome
#----------------------------------------#
my %tid2num = ();   # # of reads classified within each tax id
my $totalcount = 0;
# read tab delimited file
open(FH, "$blast") or die("Could not open $blast.\n");
foreach my $line (<FH>) {
    chomp $line;
    my ($qid, $rid, $pct, $hspl) = split("\t", $line);
    if ($pct < 95 || $hspl < 45) { next;}
    #$rid =~ /^(\S+)\_\d+\_\d+$/;
    $rid =~ /^(\S+)\_\S+\_\S+$/;
    
    $rid = $1;
    my $tid = $seq2tid{$rid};
    $tid2num{$tid}++;
    $totalcount++;
}
close FH;
#----------------------------------------#


#----------------------------------------#
# based on the abundance, pick reference
#----------------------------------------#
my $total = 0;
my %sps = ();
foreach my $tid (sort {$tid2num{$b} <=> $tid2num{$a}} keys %tid2num) {

    my $num = $tid2num{$tid}*$part;
    
    if ($num < 200) { last;}    # if abundance < 50 then ignore
    if ($total > 1000) { last;} # if more than 120 genomes have been used then stop
    
    my $sp = $tid2sp{$tid};
    if (exists $sps{$sp}) {
	if ($sps{$sp} >= 10) { next;}  # for a single species, allow max 10 reference genomes
    }
    
    $sps{$sp}++;
    $total++;

    #put the same counts for each molecule (like plasmids) of each species
    if ((($num/$totalcount)*100) >=1.0 ){
       foreach my $seq (keys %{$tid2seqs{$tid}}) {
   	  print "$seq\t$tid\t$num\t$tid2name{$tid}\n";
       }
    }
}



exit;


sub Usage {
    die("
Usage:
       perl pickrefids.pl <blast file> <partition>

Options:
       <blast>        Query sequences in FASTA format to be classified.
       <partition>    Suppose you used a quarter of total data to estimate the
                      taxonomic profile. Then partition number is 4.

Output:
       A list of NCBI genome IDs, which will be used as reference genomes.

");
}

# todo
# pick reference genomes based on the estimation of depth of coverage
# allow users to tune this parameter, right now it's hard coded
