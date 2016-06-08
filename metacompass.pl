#!/usr/bin/perl

#############################################
#
# Program: MetaCompass: Metagenomic comparative assembly
#
# Author: Bo Liu, boliu@umiacs.umd.edu
#
# Thu Aug  9 12:25:58 EDT 2012
#
#############################################

use strict;
use warnings;
use FindBin qw($Bin);


#----------------------------------------#
# read command line options
#----------------------------------------#
my $fastafile = "";
my $prefix = 0;
my $program = "";
my $niter = 3;
my $nthreads = 0;
if (scalar @ARGV == 5) {

    ($fastafile, $prefix, $program, $niter, $nthreads) = @ARGV;

} else {
    Usage();
}
#----------------------------------------#


#----------------------------------------#
# pick reference genomes
#----------------------------------------#
my $cmd = "perl $Bin/bin/pickrefseqs.pl $fastafile $prefix $nthreads";
print STDERR "$cmd\n";
system("$cmd");
my $reffile = "$prefix.refseq.fna";
#----------------------------------------#


#----------------------------------------#
# comparative assembly
#----------------------------------------#
$cmd = "perl $Bin/bin/compass.pl $fastafile $reffile $prefix $program $niter $nthreads";
print STDERR "$cmd\n";
system("$cmd");


exit;


sub Usage {
    die("
Usage:
       perl metacompass.pl <FASTA> <prefix> <mapping> <# iterations> <# threads>

Options:
       <FASTA>        DNA reads in FASTA format.
       <prefix>       Output file prefix.
       <mapping>      mummer-map or the directory of bowtie2.
       <# iterations>  # iterations to run comparative assembly. 3 is recommended.
       <# threads>    # threads to run read mapping.

Output:
      Assembled contigs:
      		prefix.[iteration #].contig
      
      Taxonomic profiling output from MetaPhyler:
      		prefix.classification; prefix.[species|genus|family|order|class|phylum].taxprof
      
      Reference genomes:
      		prefix.refseq.fna

Contact:
      Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
