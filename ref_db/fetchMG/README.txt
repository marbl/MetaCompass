
===============
 fetchMG v 1.0 
===============



fetchMG is Copyright (c) 2012 Shinichi Sunagawa, Daniel R Mende, and EMBL
(www.embl.de). fetchMG is released under the GNU General Public Licence v3.
Please see http://www.gnu.org/licenses/gpl.html and the seperately provided LICENSE file.



==============
 Introduction
==============
 
Phylogenetic markers are genes (and proteins) which can be used to reconstruct the
phylogenetic history of different organisms. One classical phylogenetic marker is the
16S ribosomal RNA gene, which is often-used but is also known to be a sub-optimal
phylogenetic marker for some organisms. Efforts to find a good set of protein coding
phylogenetic marker genes (Ciccarelli et al., Science, 2006; Sorek et al., Science, 2007)
lead to the identification of 40 universal single copy marker genes (MGs). These 40
marker genes occur in single copy in the vast majority of known organisms and they
were used to successfully reconstruct a three domain phylogenetic tree
(Ciccarelli et al., Science, 2006).



========================
 What the software does
========================
 
The program “fetchMG” was written to extract the 40 MGs from genomes and
metagenomes in an easy and accurate manner. This is done by utilizing Hidden Markov
Models (HMMs) trained on protein alignments of known members of the 40 MGs as well
as calibrated cutoffs for each of the 40 MGs. Please note that these cutoffs are only
accurate when using complete protein sequences as input files. The output of the program
are the protein sequences of the identified proteins, as well as their nucleotide sequences,
if the nucleotide sequences of all complete genes are given as an additional input.



=======
 Input
=======

A fasta file with protein coding sequences, and optionally the gene sequences of the proteins.
If the DNA sequences are available, the corresponding genes of the proteins, are also extracted.



========
 Output
========

The output of this software is saved within the specified output folder and consists of:
- 40 x COGxxxx.faa files (sequences of extracted proteins)
- 40 x COGxxxx.fna files (sequences of extracted genes)
-      marker_genes_scores.table (protein <TAB> score <TAB> marker gene ID <TAB> genome identifier)
-      temp (identifiers of proteins identified homologous to any marker gene)
-      hmmResults (specific output files from HMMer3)



=========
 Example 
=========

./fetchMG.pl -m extraction -x bin example_dataset/example_data.faa



====================
 Full program help
====================

Help & Manual
       fetchMGs.pl -h|help

Usage
       fetchMGs.pl -m|mode <extraction|calibration> [OPTIONS]

Extraction mode
       ./fetchMG.pl -m extraction <protein sequences> [optional options]
           <protein sequences>               Multi-FASTA file with protein sequences from which marker genes should be extracted
           -c|og_used                        Orthologous group id to be extracted; example: 'COG0012'; default = 'all'
           -o|outdir                         Output directory; default = 'output'
           -h|hmmdir                         Path to directory that contains hmm models; default = './lib'
           -b|bitscore                       Path to bitscore cutoff file; default = 'lib/MG_BitScoreCutoffs.defaults.txt'
           -p|protein_only                   Set if nucleotide sequences for filename.faa is not available
           -v|verybesthit_only               Recommended to use, if extracting sequences from reference genomes.
                                               For this fasta identifiers should be in the form: taxID.geneID and,
                                               if needed have ' project_id=XXX' somewhere in the header
           -t|threads                        Number of processors/threads to be used
           -d|dnaFastaFile                   Fasta file with DNA sequences of the same genes; not neccesary if protein file and dna file have the same with .faa and .fna suffixes
           -x|xbin                           Path to binaries used by this script. default = '' --> will search for variables in $PATH

Calibration mode
       ./fetchMGs.pl -m calibration <reference protein sequences> <true positives map>
           <reference protein sequences>     Multi-FASTA file with protein sequences that include marker genes (true positives)
           <true positives map>              Tab-delimited file with true positive protein identifiers and COG IDs



=======================
 Software dependencies
=======================

The fetchMG script requires the cdbyank, cdbfasta and HMMer3 executables.
These software are (c) respecitve authors, and have been installed
in the bin folder, within the fetchMG folder.

cdbyank, cdbfasta: http://sourceforge.net/projects/cdbfasta/
HMMer3:            http://hmmer.janelia.org/



=========
 Contact
=========

If you have any questions, please contact Shinichi Sunagawa.
Web: http://intranet.embl.de/research/scb/bork/members/index.php?s_personId=CP-60011932

