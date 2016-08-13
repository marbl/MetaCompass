#!/usr/bin/perl

#############################################
#
# Program: Classify sequences.
#
# Author: Bo Liu
#
# Fri Mar  9 23:24:30 EST 2012
#
#############################################

use strict;
use warnings;
use FindBin qw($Bin);

#----------------------------------------#
# read command line options
#----------------------------------------#
my $query = "";
my $blast = "blastn";
my $outdir="";
my $prefix = "mc";
my $nump = 0;
if (scalar @ARGV == 3) {
    ($query, $outdir, $nump) = @ARGV;
} else {
    Usage();
}
#----------------------------------------#


my $ref = "$Bin/markers/markers.refseq.dna";
my $param = "-word_size 28";

# run blast
my $cmd = "$blast $param -num_threads $nump -evalue 1e-10 -perc_identity 95 -outfmt 6 -max_target_seqs 100 -query $query -db $ref > $outdir/$prefix.$blast.all";
print STDERR "$cmd\n";
system("$cmd");

# filter blast
my $cmd = "$Bin/bin/best-hits.py $outdir/$prefix.$blast.all $outdir/$prefix.$blast";
print STDERR "$cmd\n";
system("$cmd");

# classification
$cmd = "$Bin/bin/metaphylerClassify $Bin/markers/markers.$blast.classifier $Bin/markers/markers.taxonomy $outdir/$prefix.$blast > $outdir/$prefix.classification";
print STDERR "$cmd\n";
system("$cmd");

$cmd = "$Bin/bin/taxprof 0.9 $outdir/$prefix.classification $outdir/$prefix $Bin/markers/tid2name.tab";
print STDERR "$cmd\n";
system("$cmd");

exit;


sub Usage {
    die("
Usage:
       perl runMetaphyler.pl <query> <prefix> <# threads>

Options:
       <query>        Query sequences in FASTA format to be classified.
       <prefix>       Output prefix.
       <# threads>    Number of threads to run BLAST.

Output:
       $prefix.blast[n/x]
                      Raw blast output.

       $prefix.classification
                      Classification results.

       $prefix.<genus|family|order|class|phylum>.taxprof
                      Taxonomy profiles at each level.

Contact:
        Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
