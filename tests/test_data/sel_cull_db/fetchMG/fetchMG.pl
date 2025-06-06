#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use List::Util qw[min max];
use FindBin;
$|++;

my $usage = "
       fetchMGs extracts the 40 single copy universal marker genes (decribed 
       in Ciccarelli et al., Science, 2006 and Sorek et al., Science, 2007)
       from genomes and metagenomes in an easy and accurate manner.

Help & Manual
       fetchMGs.pl -h|help

Usage
       fetchMGs.pl -m|mode <extraction|calibration> [OPTIONS]

Extraction mode
       ./fetchMGs.pl [options] -m extraction <protein sequences> 
           <protein sequences>           Multi-FASTA file with protein sequences from which universal single-copy marker genes should be extracted
           -o|outdir                     Output directory; default = \'output\'
           -b|bitscore                   Path to bitscore cutoff file; 
                                           default = \'\$pathInWhichThisScriptResides/lib/MG_BitScoreCutoffs.[allhits|verybesthit].txt\' (depending on -v option)
           -l|library                    Path to directory that contains hmm models; 
                                           default = \'\$pathInWhichThisScriptResides/lib\'
           -p|protein_only               Set if nucleotide sequences file for <protein sequences> is not available
           -d|dnaFastaFile               Multi-FASTA file with nucleotide sequences file for <protein sequences>; 
		                                   not neccesary if protein and nucleotide fasta file have the same name except .faa and .fna suffixes
           -v|verybesthit_only           Only extract the best hit to each OG from each genome. 
		                                   Recommended to use, if extracting sequences from reference genomes. 
                                           Please do not use for metagenomes.
                                           If this option is set fasta identifiers should be in the form: taxID.geneID and,
                                           if needed have \' project_id=XXX\' somewhere in the header
           -c|og_used                    Orthologous group id to be extracted; example: \'COG0012\'; default = \'all\'
           -t|threads                    Number of processors/threads to be used
           -x|executables                Path to executables used by this script (hmmsearch; cdbfasta, cdbyank). 
		                                   default = \'\$pathInWhichThisScriptResides/bin\' 
		                                   If set to \'\' will search for executables in \$PATH

Calibration mode
       ./fetchMGs.pl [options] -m calibration <reference protein sequences> <true positives map>
           <reference protein sequences> Multi-FASTA file with protein sequences that include marker genes (true positives)
           <true positives map>          Tab-delimited file with true positive protein identifiers and COG IDs

           -o|outdir                     Output directory; default = \'calibration\'		   
           -b|bitscore                   Path to bitscore cutoff file; 
                                           default = \'\$pathInWhichThisScriptResides/lib/MG_BitScoreCutoffs.uncalibrated.txt\'
           -v|verybesthit_only           Use if calibrating for extraction using the -v option 
										   
           The other options for \'-m extraction\' can also be used here.
";



################################
#Defaults / Global variables
my $bin                     = "DEFAULT";
my $hmm_dir                 = "DEFAULT";
my $mode                    = "";
my $cog_used                = "all";
my $outdir                  = "DEFAULT";
my $cutoff_file             = "DEFAULT";
my $protein_only            = 0;
my $besthit_only            = 0;
my $intProcessorNumber      = 1;
my $floatMin                = 60;
my $manual                  = 0;

my $strProteinFastaFileName = "";
my $strDNAFastaFileName     = "";
my $strPrefix               = "test";
my $mapFile;

my @arrAllResults;
my %hashCutoffs;

#replace the contents of this array if less/different COGs should be used
my @arrCOGsIDs = ( 'COG0012', 'COG0016', 'COG0018', 'COG0048', 'COG0049', 'COG0052', 'COG0080', 'COG0081', 'COG0085', 'COG0087', 'COG0088', 'COG0090', 'COG0091', 'COG0092', 'COG0093', 'COG0094', 'COG0096', 'COG0097', 'COG0098', 'COG0099', 'COG0100', 'COG0102', 'COG0103', 'COG0124', 'COG0172', 'COG0184', 'COG0185', 'COG0186', 'COG0197', 'COG0200', 'COG0201', 'COG0202', 'COG0215', 'COG0256', 'COG0495', 'COG0522', 'COG0525', 'COG0533', 'COG0541', 'COG0552' );

################################
#Get options

GetOptions(
	'm|mode=s'				=> \$mode,
	'h|help'				=> \$manual,
	'c|og_used=s'			=> \$cog_used,
	'o|outdir=s'			=> \$outdir,
	'l|library=s'			=> \$hmm_dir,
	'b|bitscore=s'			=> \$cutoff_file,
	'p|protein_only'		=> \$protein_only,
	'v|verybesthit_only'	=> \$besthit_only,
	't|threads=i'			=> \$intProcessorNumber,
	'd|dnaFastaFile=s'		=> \$strDNAFastaFileName,
	'x|executables=s'		=> \$bin
);
################################

if ( $mode eq 'calibration' ) {
	if ( $outdir eq 'DEFAULT' ) {
		$outdir = "calibration";
	}
}
elsif ( $mode eq 'extraction' ) {
	if ( $outdir eq 'DEFAULT' ) {
		$outdir = "output";
	}
}

if ( $mode eq 'calibration' ) {
	if ( $cutoff_file eq 'DEFAULT' ) {
		$cutoff_file = "$FindBin::Bin/lib/MG_BitScoreCutoffs.uncalibrated.txt";
	}
}
elsif ( $mode eq 'extraction' ) {
	if ($besthit_only) {
		if ( $cutoff_file eq 'DEFAULT' ) {
			$cutoff_file = "$FindBin::Bin/lib/MG_BitScoreCutoffs.verybesthit.txt";
		}
	}
	else {
		if ( $cutoff_file eq 'DEFAULT' ) {
			$cutoff_file = "$FindBin::Bin/lib/MG_BitScoreCutoffs.allhits.txt";
		}
	}
}

if ( $bin eq 'DEFAULT' ) {
	$bin = "$FindBin::Bin/bin/";
}

if ( $hmm_dir eq 'DEFAULT' ) {
	$hmm_dir = "$FindBin::Bin/lib/";
}


#check bin directory format
if ( $bin ne "") {
	if ( not $bin =~ m/\/$/ ) {
		$bin = $bin . "/";
	}
}

#check fasta file input and extract prefix as run name
if ($manual) {
	pod2usage( -exitstatus => 0, -verbose => 2 );
	exit 0;
}

if ( ( $mode eq "extraction" ) || ( $mode eq "calibration" ) ) {
	if ( ( not ( $ARGV[0] ) ) or ( not ( -e $ARGV[0] ) ) ) {
		die $usage . "ERROR: Missing <input protein fasta file>\n";
	}
	else {
		$strProteinFastaFileName = $ARGV[0];
		
		#protein file ends in .faa --> automatically generate $strDNAFastaFileName if not set before
		if ( $strProteinFastaFileName =~ /.faa$/ ) {
			my @file_array = split( "\/", $strProteinFastaFileName );

			my $filenameWithoutPath = pop(@file_array);
			$filenameWithoutPath =~ m/(.+)\.(faa)/;
			$strPrefix = $1;

			if ( length($strDNAFastaFileName) == 0 ) {
				if ( ( scalar @file_array ) > 0 ) {
					my $strPathToFastaFiles = join( "\/", @file_array );
					$strDNAFastaFileName = $strPathToFastaFiles . "\/" . $strPrefix . ".fna";
				}
				else {
					$strDNAFastaFileName = $strPrefix . ".fna";
				}
			}
		}

		#protein file does not end in .faa --> need $strDNAFastaFileName from option -d
		else {
			if ( length($strDNAFastaFileName) == 0 && $protein_only == 0 ) {
				die $usage . "ERROR: Missing input dna fasta file; please specify it using -d\n";
			}

			my @file_array = split( "\/", $strProteinFastaFileName );

			my $filenameWithoutPath = pop(@file_array);

			if ( $filenameWithoutPath =~ /\./ ) {
				my @file_arrayWOPath = split( /\./, $filenameWithoutPath );
				pop(@file_arrayWOPath);
				$strPrefix = join( "\.", @file_arrayWOPath );
			}
			else {
				$strPrefix = $filenameWithoutPath;
			}
		}
	}
	# check if a file named $strDNAFastaFileName exist; if not die
	
	if ( ( $protein_only == 0 ) &&  ( not (-e $strDNAFastaFileName ))) {
		die $usage . "ERROR: Missing input dna fasta file\n";
	}
	
}


print "
===================================================================================
         FetchMGs v1.0 - extraction of marker genes from protein sequences
          Copyright (c) 2012 Shinichi Sunagawa, Daniel R Mende, and EMBL
===================================================================================

";

if ( $mode eq "extraction" ) {
	&extractMGs();
}
elsif ( $mode eq "calibration" ) {

	$mapFile = $ARGV[1];
	&extractMGs();

	print "\nCalibrating cutoffs...";
	&calibrateBitScores( $strPrefix, $mapFile, $outdir );
	if ($besthit_only) {
		print "done\n\nTo use the new calibrations, run fetchMGs with option: -b $outdir\/MG_BitScoreCutoffs.verybesthit.txt\n";
	}
	else {
		print "done\n\nTo use the new calibrations, run fetchMGs with option: -b $outdir\/MG_BitScoreCutoffs.txt\n";
	}
}
else {
	die $usage;
}

#SUBROUTINES
# global variables needed:
# $strProteinFastaFileName
# $strDNAFastaFileName
# $cutoff_file
sub extractMGs {

	if ( $cog_used ne "all" ) {
		@arrCOGsIDs = ( $cog_used );
	}

	#Prepare output directory
	if ( -e $outdir ) {
		system "rm -r $outdir;";
		mkdir $outdir;
	}
	else {
		mkdir $outdir;
	}
	print "Parsing cutoffs... ";
	&parseCutoffsFile($cutoff_file);
	print "done\n";

	print "Running HMMsearch... ";
	&allCOGsHMMsearch($strProteinFastaFileName);
	print "done\n";

	print "Processing results...\n";
	&filterResultsArray_besthitAmongCOGs();

	if ($besthit_only) {
		&filterResultsArray_besthitAmongGenomes();
	}

	&printSequences( $strProteinFastaFileName, $strDNAFastaFileName, $strPrefix );
	print "Processing completed!\n\nExtracted sequences are saved in the directory: $outdir\n";
}

#SUBROUTINES
# global variables needed:
# %hashCutoffs
# $besthit_only
# $mode
sub parseCutoffsFile {
	my $cutoff_file = $_[0];

	open( CUTOFFS, "<$cutoff_file" ) or die "ERROR: Cannot open $cutoff_file";

	chomp( my $header = <CUTOFFS> );
	unless ( $header =~ m/^#/ ) {
		die "\nERROR: $cutoff_file is not a valid cutoff file for fetchMGs";
	}

	if ( $mode eq "extraction" ) {
		if ($besthit_only) {
			unless ( $header eq "#CALIBRATED CUTOFFS FILE - BEST HITS" ) {
				die "\nERROR: To extract using the option v|verybesthit, the file has to be calibrated with the -v option";
			}
		}
		else {
			unless ( $header eq "#CALIBRATED CUTOFFS FILE - ALL HITS" ) {
				die "\nERROR: To extract without the option v|verybesthit, the file has to be calibrated without the -v option";
			}
		}
	}
	elsif ( $mode eq "calibration" ) {
		unless ( $header eq "#UNCALIBRATED CUTOFFS FILE" ) {
			die "\nERROR: When calibrating, the default bit score file (provided with fetchMGs) has to be used";
		}
	}
	while ( my $currentCutoffLine = <CUTOFFS> ) {
		chomp($currentCutoffLine);
		my @arrCurrentCutoffLine = split( "\t", $currentCutoffLine );
		my $strCOG               = $arrCurrentCutoffLine[0];
		my $cutoff               = $arrCurrentCutoffLine[1];

		$hashCutoffs{$strCOG} = $cutoff;
	}
	close(CUTOFFS);
}

#run HMMer3 - internal subroutine
# global variables needed:
# none
sub runHmmerSearch {
	my $strCOGid      = $_[0];
	my $floatCutoff   = $_[1];
	my $strHMMfile    = $_[2];
	my $strFastaFile  = $_[3];
	my $strRunName    = $_[4];
	my $strOutfolder  = $_[5];
	my $intProcessors = $_[6];

	my $strDomainOutput = $strOutfolder . "\/" . $strRunName . "." . $strCOGid . ".dom";
	my $strHmmerOutput  = $strOutfolder . "\/" . $strRunName . "." . $strCOGid . ".out";

	my $cmd = $bin . "hmmsearch --noali --cpu " . $intProcessors . " -o " . $strHmmerOutput . " --domtblout " . $strDomainOutput . " -T " . $floatCutoff . " " . $strHMMfile . " " . $strFastaFile;

	system($cmd);
	return ($strDomainOutput);
}

#parse HMMer3 results - internal subroutine
# global variables needed:
# none
sub parseHmmerSearch {
	my $strCOGid        = $_[0];
	my $strDomainOutput = $_[1];
	my $floatCutoff     = $_[2];

	my @arrResults = ();
	open( INFILE, "<$strDomainOutput" ) or die "Cannot open the infile file for data\n";

	while ( my $currentLine = <INFILE> ) {
		chomp($currentLine);

		if ( $currentLine =~ m/^#/ ) {
			next;
		}

		my @currentLineArray = split( /\s+/, $currentLine );
		my $strQueryId       = $currentLineArray[0];
		my $floatBitscore    = $currentLineArray[7];

		my @arrQuery = split( /\./, $strQueryId );
		my $strTaxID = $arrQuery[0];

		#get project ID if encoded in fasta header / description
		my $strProjectID = "_NA_";
		for my $element (@currentLineArray) {
			next unless $element =~ m/^project_id=/;
			my @element_array = split( /=/, $element );
			$strProjectID = $element_array[1];
			$strProjectID =~ s/"//g;
		}

		my $strTaxProjectID = $strTaxID . "." . $strProjectID;

		my @tupleResultsLine = ( $strQueryId, $strTaxProjectID, $strCOGid, $floatBitscore );

		#if score higher than cutoff add to array of arrays of results (@arrResults)
		if ( $floatBitscore >= $floatCutoff ) {
			push @arrResults, [@tupleResultsLine];
		}
	}
	return ( \@arrResults );
}

#run and parse HMMer3 results
# global variables needed:
# $outdir
# @arrCOGsIDs
# $hmm_dir
# %hashCutoffs
# $intProcessorNumber
# @arrAllResults

sub allCOGsHMMsearch {
	$strProteinFastaFileName = $_[0];

	my $strHMMoutdir      = $outdir . "\/hmmResults";
	my $strCurrentRunName = "markerGenes";
	mkdir $strHMMoutdir;

	foreach my $strCurrentCOG (@arrCOGsIDs) {

		my $strCurrentHMM      = $hmm_dir . "\/" . $strCurrentCOG . ".hmm";
		my $floatCurrentCutoff = $hashCutoffs{$strCurrentCOG};

		my $strDomainOutput = &runHmmerSearch( $strCurrentCOG, $floatCurrentCutoff, $strCurrentHMM, $strProteinFastaFileName, $strCurrentRunName, $strHMMoutdir, $intProcessorNumber );
		my $refArrResults = &parseHmmerSearch( $strCurrentCOG, $strDomainOutput, $floatCurrentCutoff );

		push( @arrAllResults, @{$refArrResults} );
	}

}

#filtering of results
#1. besthit among cogs
#
# global variables needed:
# @arrAllResults
sub filterResultsArray_besthitAmongCOGs {

	my %hashGene2Results;
	my %hashGene2Keep;
	
	for my $i ( 0 .. $#arrAllResults ) {

		my $strQueryId      = $arrAllResults[$i][0];
		my $strTaxProjectID = $arrAllResults[$i][1];
		my $strCOGid        = $arrAllResults[$i][2];
		my $floatBitscore   = $arrAllResults[$i][3];

		if ( exists( $hashGene2Results{$strQueryId} ) ) {

			my $floatPrevBitScore = $hashGene2Results{$strQueryId}[3];

			if ( $floatBitscore > $floatPrevBitScore ) {
				
				$hashGene2Results{$strQueryId} = [ $strQueryId, $strTaxProjectID, $strCOGid, $floatBitscore ];
				$hashGene2Keep{$strQueryId} =  1;
			}
			elsif ( $floatBitscore == $floatPrevBitScore ) {
				if ($hashGene2Results{$strQueryId}[2] ne $strCOGid){
				
					$hashGene2Keep{$strQueryId} =  0;
				}
			}
		}
		else {
			$hashGene2Results{$strQueryId} = [ $strQueryId, $strTaxProjectID, $strCOGid, $floatBitscore ];
			$hashGene2Keep{$strQueryId} =  1;
		}
	}

	@arrAllResults = ();
	for my $refCurrentResult ( sort keys %hashGene2Results ) {
		if ($hashGene2Keep{$refCurrentResult}){
			push( @arrAllResults, [ @{ $hashGene2Results{$refCurrentResult} } ] );
		}
		#else{
		#	print "same score for 2 OGs, removing:\t$refCurrentResult\n";
		#}
	}
}

#filtering of results
# 2. besthit per genomes
#
# global variables needed:
# @arrAllResults
sub filterResultsArray_besthitAmongGenomes {

	my %hashOG2Project2Results;
	
	for my $i ( 0 .. $#arrAllResults ) {

		my $strQueryId      = $arrAllResults[$i][0];
		my $strTaxProjectID = $arrAllResults[$i][1];
		my $strCOGid        = $arrAllResults[$i][2];
		my $floatBitscore   = $arrAllResults[$i][3];

		if ( exists( $hashOG2Project2Results{$strCOGid} ) and exists( $hashOG2Project2Results{$strCOGid}{$strTaxProjectID} ) ) {

			my $floatPrevBitScore = $hashOG2Project2Results{$strCOGid}{$strTaxProjectID}[3];

			if ( $floatBitscore > $floatPrevBitScore ) {
				$hashOG2Project2Results{$strCOGid}{$strTaxProjectID} = [ $strQueryId, $strTaxProjectID, $strCOGid, $floatBitscore ];
				
			}
		}
		else {
			$hashOG2Project2Results{$strCOGid}{$strTaxProjectID} = [ $strQueryId, $strTaxProjectID, $strCOGid, $floatBitscore ];
		}
	}

	@arrAllResults = ();
	for my $strCurrentOG ( sort keys %hashOG2Project2Results ) {
		for my $strCurrentProject ( sort keys %{ $hashOG2Project2Results{$strCurrentOG} } ) {
			push( @arrAllResults, [ @{ $hashOG2Project2Results{$strCurrentOG}{$strCurrentProject} } ] );
		}
	}
}

#print sequences and make MG_SCORES files for calibration
#
# download faSomeRecords from:
# http://hgdownload.cse.ucsc.edu/admin/exe/
#
# global variables needed:
# @arrAllResults
# @arrCOGsIDs
# $outdir
# $protein_only
sub printSequences {

	my $strProteinFastaFileName = $_[0];
	my $strDNAFastaFileName     = $_[1];
	my $strPrefix               = $_[2];

	#indexing using cdbfasta
	my $cmd = $bin . "cdbfasta $strProteinFastaFileName";
	system $cmd;

	if ( $protein_only == 0 ) {
		my $cmd = $bin . "cdbfasta $strDNAFastaFileName";
		system $cmd;
	}

	#make ID and scores files
	my $strIDFilesOutdir = $outdir . "\/temp";
	mkdir $strIDFilesOutdir;
	foreach my $strCurrentCOG (@arrCOGsIDs) {
		open( ID_OUT, ">$strIDFilesOutdir/$strCurrentCOG.IDs.txt" );
		close(ID_OUT);
	}

	my $marker_genes_scores_table = "$outdir\/$strPrefix.$cog_used.marker_genes_scores.table";
	open( MG_SCORES, ">$marker_genes_scores_table" );

	print MG_SCORES "#protein_sequence_id\tHMM bit score\tCOG\ttaxid.projectid (if possible)\n";

	for my $i ( 0 .. $#arrAllResults ) {

		my $strQueryId      = $arrAllResults[$i][0];
		my $strTaxProjectID = $arrAllResults[$i][1];
		my $strCurrentCOG   = $arrAllResults[$i][2];
		my $floatBitscore   = $arrAllResults[$i][3];

		open( ID_OUT, ">>$strIDFilesOutdir/$strCurrentCOG.IDs.txt" );
		print ID_OUT $strQueryId . "\n";
		close(ID_OUT);

		my $mg_scores_line = $strQueryId . "\t" . $floatBitscore . "\t" . $strCurrentCOG . "\t" . $strTaxProjectID . "\n";
		print MG_SCORES $mg_scores_line;
	}
	close(MG_SCORES);

	foreach my $strCurrentCOG (@arrCOGsIDs) {

		my $cmd = "cat $strIDFilesOutdir/$strCurrentCOG.IDs.txt | " . $bin . "cdbyank $strProteinFastaFileName.cidx > ${outdir}\/${strCurrentCOG}.faa";

		system $cmd;

		if ( $protein_only  == 0 ) {
			$cmd = "cat $strIDFilesOutdir/$strCurrentCOG.IDs.txt | " . $bin . "cdbyank $strDNAFastaFileName.cidx > ${outdir}\/${strCurrentCOG}.fna";

			system $cmd;
		}
	}
}

#####
#Calibrate bitscores using results from the extraction and a mapping file of known members of OGs in the fastafile
# global variables needed:
# @arrCOGsIDs
# $outdir
# $floatMin
# $cog_used
sub calibrateBitScores {
#################################################
	#Get maximum F-scores for marker gene predictions
	#Precision = p = (tp / (tp + fp))
	#Recall    = r = (tp / (tp + fn))
	#F-score   = (1 + beta^2) (p*r) / ((beta^2 * p) +r))
	#$beta = 1;

	my $strPrefix      = $_[0];
	my $strMappingFile = $_[1];

	#initialize variables
	my $beta = 1;
	my %hashGene2COG;
	my %hashCOG2Count;
	my %hashCOG2PositiveBitScoreArray;
	my %hashCOG2NegativeBitScoreArray;

	foreach my $strCOG (@arrCOGsIDs) {
		$hashCOG2Count{$strCOG}                 = 0;
		$hashCOG2PositiveBitScoreArray{$strCOG} = [];
		$hashCOG2NegativeBitScoreArray{$strCOG} = [];
	}

	#read in mapping files gene --> OG (expected to be positive/ training data)
	open( MAP, "<$strMappingFile" );
	while ( my $currentMapLine = <MAP> ) {
		chomp($currentMapLine);
		my @arrCurrentMapLine = split( "\t", $currentMapLine );
		my $geneName          = $arrCurrentMapLine[0];
		my $strCOG            = $arrCurrentMapLine[1];

		$hashGene2COG{$geneName} = $strCOG;
		$hashCOG2Count{$strCOG}++;
	}
	close(MAP);

	#read in scores tables
	my $marker_genes_scores_table = "$outdir\/$strPrefix.$cog_used.marker_genes_scores.table";
	print("\nProcessing $marker_genes_scores_table \n");
	open( MG_SCORES, "<$marker_genes_scores_table" );

	while ( my $currentScoreLine = <MG_SCORES> ) {
		if ( !( $currentScoreLine =~ m/^#/ ) ) {
			chomp($currentScoreLine);
			my @arrCurrentScoreLine = split( "\t", $currentScoreLine );

			my $geneName = $arrCurrentScoreLine[0];
			my $bitScore = $arrCurrentScoreLine[1];
			my $strCOG   = $arrCurrentScoreLine[2];

			if ( exists( $hashGene2COG{$geneName} ) and $strCOG eq $hashGene2COG{$geneName} ) {
				push( @{ $hashCOG2PositiveBitScoreArray{$strCOG} }, $bitScore );
			}
			else {
				push( @{ $hashCOG2NegativeBitScoreArray{$strCOG} }, $bitScore );
			}
		}
	}
	close(MAP);

	#calculate F-Scores for all used OGs
	print("Generating new cutoffs file... ");

	if ($besthit_only) {
		open( OUT, ">$outdir\/MG_BitScoreCutoffs.verybesthit.txt" );
	}
	else {
		open( OUT, ">$outdir\/MG_BitScoreCutoffs.txt" );
	}

	if ($besthit_only) {
		print OUT "#CALIBRATED CUTOFFS FILE - BEST HITS\n";
	}
	else {
		print OUT "#CALIBRATED CUTOFFS FILE - ALL HITS\n";
	}

	foreach my $strCOG (@arrCOGsIDs) {

		open( OUT2, ">$outdir/$strCOG.bitscore.result" );

		my @arrPositiveBitScores = @{ $hashCOG2PositiveBitScoreArray{$strCOG} };
		my @arrNegativeBitScores = @{ $hashCOG2NegativeBitScoreArray{$strCOG} };

		my @arrPositiveBitScoresSorted = sort { $a <=> $b } @arrPositiveBitScores;
		my @arrNegativeBitScoresSorted = sort { $a <=> $b } @arrNegativeBitScores;

		my $max_Fscore_bitscore = 0;
		my $max_Fscore          = 0;
		my $max_Precision       = 0;
		my $max_TruePos         = 0;
		my $max_FalsePos        = 0;
		my $max_FalseNeg        = 0;
		my $max_Recall          = 0;

		my $intMin = int( $floatMin / 10.0 ) * 10.0;
		my $intMax = 0;

		if ( scalar @arrPositiveBitScoresSorted == 0 and scalar @arrNegativeBitScoresSorted == 0 ) {
			$intMax = $intMin;
		}
		elsif ( scalar @arrPositiveBitScoresSorted == 0 ) {
			$intMax = $arrNegativeBitScoresSorted[-1];
		}
		elsif ( scalar @arrNegativeBitScoresSorted == 0 ) {
			$intMax = $arrPositiveBitScoresSorted[-1];
		}
		else {
			$intMax = int( max( $arrPositiveBitScoresSorted[-1], $arrNegativeBitScoresSorted[-1] ) ) + 1;
		}

		$max_Fscore_bitscore = $intMin;

		# go through all possible cutoffs
		for ( my $intBitCutoff = $intMin ; $intBitCutoff <= $intMax ; $intBitCutoff += 10 ) {
			my $intTruePos  = 0;
			my $intFalsePos = 0;
			my $intFalseNeg = 0;

			#throw away all values below cutoff
			while ( scalar(@arrPositiveBitScoresSorted) > 0 and $arrPositiveBitScoresSorted[0] < $intBitCutoff ) {
				shift(@arrPositiveBitScoresSorted);
			}

			while ( scalar(@arrNegativeBitScoresSorted) > 0 and $arrNegativeBitScoresSorted[0] < $intBitCutoff ) {
				shift(@arrNegativeBitScoresSorted);
			}

			$intTruePos  = scalar(@arrPositiveBitScoresSorted);
			$intFalseNeg = $hashCOG2Count{$strCOG} - $intTruePos;
			$intFalsePos = scalar(@arrNegativeBitScoresSorted);
			my $floatPrecision = 0.0;
			my $floatRecall    = 0.0;
			my $floatFscore    = 0.0;
			if ( $intTruePos != 0 ) {
				$floatPrecision = $intTruePos / ( $intTruePos + $intFalsePos );
				$floatRecall    = $intTruePos / ( $intTruePos + $intFalseNeg );
				$floatFscore = ( 1 + ( $beta * $beta ) ) * ( $floatPrecision * $floatRecall ) / ( ( $beta * $beta * $floatPrecision ) + $floatRecall );
			}

			if ( $floatFscore >= $max_Fscore ) {
				$max_Fscore_bitscore = $intBitCutoff;
				$max_Fscore          = $floatFscore;
				$max_Precision       = $floatPrecision;
				$max_TruePos         = $intTruePos;
				$max_FalsePos        = $intFalsePos;
				$max_FalseNeg        = $intFalseNeg;
				$max_Recall          = $floatRecall;
			}

			print OUT2 "Bitscore: $intBitCutoff\tDetected: " . ( $intTruePos + $intFalsePos ) . "\tTrue: $intTruePos\tFalse: $intFalsePos\tMissing: $intFalseNeg\tPrecision: $floatPrecision\tRecall: $floatRecall\tFscore: $floatFscore\n";
		}

		print OUT2 "Bitscore for max. score: $max_Fscore_bitscore\tDetected: " . ( $max_TruePos + $max_FalsePos ) . "\tTrue: $max_TruePos\tFalse: $max_FalsePos\tMissing: $max_FalseNeg\tPrecision: $max_Precision\tRecall: $max_Recall\tFscore: $max_Fscore\n";
		close OUT2;

		print OUT "$strCOG\t$max_Fscore_bitscore\t$max_Precision\t$max_Recall\n";
	}
	close OUT;
}

__END__


=head1 Name

B<fetchMGs>

=head1 Version

1.0

=head1 Software description

fetchMGs extracts the 40 single copy universal marker genes (decribed in 
Ciccarelli et al., Science, 2006 and Sorek et al., Science, 2007) from genomes and 
metagenomes in an easy and accurate manner. This is done by utilizing profile Hidden Markov
Models (HMMs) trained on protein alignments of known members of the 40 MGs as well
as calibrated cutoffs for each of the 40 MGs. Please note that these cutoffs are only
accurate when using complete protein sequences as input files. The output of the program
are the protein sequences of the identified proteins, as well as their nucleotide sequences,
if the nucleotide sequences of all complete genes are given as an additional input.

=head1 Input

A fasta file with protein coding sequences, and optionally the nucleotide sequences of the proteins. If the DNA sequences are available, the corresponding nucleotide sequences of the proteins, are also extracted.

=head1 Output

The output of this software is saved within the specified output folder and consists of:

=over

=item - 40 x COGxxxx.faa files (sequences of extracted proteins)

=item - 40 x COGxxxx.fna files (sequences of extracted genes)

=item      - marker_genes_scores.table (protein <TAB> score <TAB> marker gene ID <TAB> genome identifier)

=item      - temp (identifiers of proteins identified homologous to any marker gene)

=item      - hmmResults (specific output files from HMMer3)

=back

=head1 Synopsis

fetchMGs.pl -m|mode <I<extraction>|I<calibration>> [OPTIONS]



=head1 Extraction mode

./fetchMGs.pl [options] -m extraction <protein sequences> 

I<Required options>

=over

=over

=item <protein sequences>

Multi-FASTA file with protein sequences from which marker genes should be extracted

=back

=back

I<Further options>


=over

=over

=item -o|outdir

Output directory; default = "output"

=back

=back


=over

=over

=item -b|bitscore

Path to bitscore cutoff file; Path to bitscore cutoff file; 
default = "$pathInWhichThisScriptResides/lib/MG_BitScoreCutoffs.[allhits|verybesthit].txt" (depending on -v option)

=back

=back


=over

=over

=item -l|library

Path to directory that contains hmm models; 
default = "$pathInWhichThisScriptResides/lib"

=back

=back


=over

=over

=item -p|protein_only

Set if nucleotide sequences file for <protein sequences> is not available

=back

=back


=over

=over

=item -d|dnaFastaFile

Multi-FASTA file with nucleotide sequences file for <protein sequences>; 
Not neccesary if protein and nucleotide fasta file have the same name except .faa and .fna suffixes

=back

=back


=over

=over

=item -v|verybesthit_only 

Only extract the best hit to each OG from each genome.
Recommended to use, if extracting sequences from reference genomes. 
Please do not use for metagenomes.
If this option is set fasta identifiers should be in the form: taxID.geneID and, if needed, have " project_id=XXX" somewhere in the header

=back

=back


=over

=over

=item -c|og_used

Orthologous group id to be extracted; example: "COG0012"; default = "all"

=back

=back


=over

=over

=item -t|threads

Number of processors/threads to be used

=back

=back


=over

=over

=item -x|executables

Path to binaries used by this script. default = "" --> will search for variables in \$PATH
Path to executables used by this script (hmmsearch; cdbfasta, cdbyank). 
default = "$pathInWhichThisScriptResides/bin" 
If set to "" will search for executables in \$PATH
										   
=back

=back



=head1 Calibration mode

./fetchMGs.pl -m calibration <reference protein sequences> <true positives map>

I<Required options>

=over

=over

=item <reference protein sequences>

Multi-FASTA file with protein sequences that include marker genes (true positives)

=back

=back

I<Further options>

=over

=over

=item <true positives map>

Tab-delimited file with true positive protein identifiers and COG IDs

=back

=back


=over

=over

=item -o|outdir

Output directory; default = "output"

=back

=back


=over

=over

=item -b|bitscore

Path to bitscore cutoff file; Path to bitscore cutoff file; 
default = "$pathInWhichThisScriptResides/lib/MG_BitScoreCutoffs.uncalibrated.txt"

=back

=back

=over

=over

=item -v|verybesthit

Use if calibrating for extraction using the -v option

=back

=back


=over

The other options for \'-m extraction\' can also be used here.

=back



=head1 Software dependencies

The fetchMGs script requires the cdbyank, cdbfasta and HMMer3 executables.
These software are (c) respecitve authors, and have been installed
in the bin folder, within the fetchMGs folder. If these are incompatible
with your system please install them yourself.

cdbyank, cdbfasta:  http://sourceforge.net/projects/cdbfasta/

hmmsearch (HMMer3): http://hmmer.janelia.org/


=head1 Author

The B<fetchMGs> package was developed by Shinichi Sunagawa and Daniel R Mende (Bork Group, EMBL) (http://www.bork.embl.de). External software used by the B<fetchMGs> package are copyright respective authors. 

=head1 Copyright

Copyright (c) 2012 Shinichi Sunagawa, Daniel R Mende, and EMBL. fetchMGs is released under the GNU General Public Licence v3 (http://www.gnu.org/licenses/gpl.html). 


