#!/usr/bin/perl

#===============================================================================
# Author: Victoria Cepeda (vcepeda@umiacs.umd.edu)
#
# File: metacompass.pl
# Date: 2015-09-28
# Version:1.0
#
# Description: Metagenomic comparative assembly tool
#
#===============================================================================

use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;

#use File::Copy qw(move);

# read command line options
my (
    $fastafile, $fastqfile, $reffile, $outdir, $niter,
    $nthreads,  $pickref,   $mincov,  $minlen, $help
);

GetOptions(
    'f|fasta:s'         => \$fastafile,
    'q|fastq:s'         => \$fastqfile,
    'r|refseq:s'        => \$reffile,
    'o|output_dir:s'    => \$outdir,
    'ni|n_iterations:i' => \$niter,
    'nt|n_threads:i'    => \$nthreads,
    'c|mincov:i'        => \$mincov,
    'l|minlen:i'        => \$minlen,
    'k|pickref:s'       => \$pickref,
    'h|help:s'          => \$help,
);
Usage() unless ( $fastafile || $fastqfile );
$outdir   ||= "MetaCompass_out";
$niter    ||= 1;
$nthreads ||= 1;
$pickref  ||= "breadth";
$mincov   ||= "2";
$minlen   ||= "300";

mkdir($outdir) unless ( -d $outdir );
my $prefix = "mc";

#create log file
my $filelog = "$outdir/log.txt";
open( my $outfile, '>', $filelog );
print $outfile "MetaCompassv1.0\n";

#Convert input to fasta
if ( !$fastafile ) {
    $fastafile = "$outdir/reads.fasta";
    print $outfile "#Converting fastq to fasta \n";
    my $cmd = "perl $Bin/bin/fq2fa.pl -i $fastqfile -o $fastafile";
    print $outfile "$cmd\n";
    system("$cmd   1>>$filelog 2>&1");
}

#Select references
if ( !$reffile ) {
    &pick_reference_genomes( $fastafile, $outdir, $nthreads );
    $reffile = "$outdir/mc.refseq.fna";
}

#Assembly
&comparative_assembly(
    $fastafile, $reffile, $outdir, $niter,
    $nthreads,  $pickref, $mincov, $minlen
);

# Select reference genomes using MetaPhyler
sub pick_reference_genomes {
    print $outfile "##Selecting reference genomes\n";
    my $cmd = "perl $Bin/bin/pickrefseqs.pl $_[0] $_[1] $_[2]";
    print $outfile "$cmd\n";
    system("$cmd   1>>$filelog 2>&1");
}

#Build contigs using comparative assembly
sub comparative_assembly {
    print $outfile "##Build contigs using comparative assembly\n";
    for ( my $i = 1 ; $i <= $niter ; ++$i ) {

        print $outfile "### iteration ", $i, " ###\n";

        if ( $i != 1 ) {
            my $j = $i - 1;
            $reffile = "$outdir/$prefix.newref_$j.fasta";
        }

        my $cmd = "bowtie2-build -q $reffile $outdir/$prefix.refseq.$i";
        print $outfile "# Build bowtie2 index\n";
        print $outfile "$cmd\n";
        system("$cmd   1>>$filelog 2>&1");

        $cmd =
"bowtie2 --sam-nohead --sam-nosq --end-to-end --quiet -k 30 -p $nthreads -x $outdir/$prefix.refseq.$i -f $fastafile > $outdir/$prefix.$i.bowtie2";
        print $outfile "# Run bowtie2 read mapping\n";
        print $outfile "$cmd\n";
        system("$cmd");

        if ( $i != 1 ) { $pickref = "all"; }
        $cmd =
"$Bin/bin/buildcontig -s $outdir/$prefix.$i.bowtie2 -r $reffile -o $outdir -c $mincov -l $minlen -n T -k $pickref";
        print $outfile "# Build contigs\n";
        print $outfile "$cmd\n";
        system("$cmd   1>>$filelog 2>&1");
        rename "$outdir/contigs.fasta", "$outdir/contigs_$i.fasta";
        rename "$outdir/newref.fasta",  "$outdir/mc.newref_$i.fasta";
        print "\n";
    }

}
exit;

sub Usage {
    die( "
Program: MetaCompass (Metagenomic comparative assembly tool)
Version: 1.0 

Usage:
        perl metacompass.pl -f/-q <FASTA/FASTQ> [options]

        <FASTA/FASTQ>	DNA reads in FASTA or FASTQ format.

Options:
	-o/--output <dir>		Output directory (default: MetaCompass_out).
	-ni/--num_iterations <int>	#iterations to run comparative assembly. 3 is recommended (default: 1)
       	-nt/--num_threads <int>		#threads to run metaphyler and read mapping (default: 1).
	-r/--refseq <FASTA>		reference sequences used to guide genome assembly
        -c/--mincov <int>		minimum depth of coverage for contigs (default: 2)
        -l/--minlen <int>		minimum length for contigs (default: 300bp)
        -k/--pickref <string>		pick reference if there are multiple best ones (default: breadth)
		      all	use all of them; can produce redundant contigs
                      breadth	pick reference with highest breadth of coverage
                      depth	pick reference with highest depth of coverage

" );
}
