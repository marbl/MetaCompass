# MetaCompassv1.0
Last updated: March 30th, 2016



# REQUIREMENTS:
GNU C/C++; Perl; BLAST; Bowtie 2
#todo:add software


# INSTALLATION:
    git clone https://github.com/marbl/MetaCompass.git
    cd MetaCompass
    ./install.pl


# USAGE

-- I have a set of metagenomic reads, and want to perform comparative assembly.

    ./metacompass.pl -f/-q <input FASTA/FASTQ reads> [options] 


-- I know the reference genomes, or I want to perform comparative assembly for a particular genome.

    ./metacompass.pl -f/-q [input FASTA/FASTQ reads] -r [reference genomes in FASTA] [options]


-- You can try MetaCompass on two test data sets.
For a small data set assuming reference genomes are known:

    ./metacompass.pl -f test/small.reads.fna -r test/small.refgeno.fna -o test1_out

For a big data set assuming reference genomes are known:

    ./metacompass.pl -f test/big.reads.fna -r test/big.refgeno.fna -o test2_out

For a big data set assuming reference genomes are NOT known:

    ./metacompass.pl -f test/big.reads.fna -o test3_out

# Output:
      Assembled contigs:
                contigs_[iteration #].fasta
      
      Taxonomic profiling output from MetaPhyler:
                mc.classification; mc.[species|genus|family|order|class|phylum].taxprof
      
      Selected Reference genomes:
                mc.refseq.fasta
