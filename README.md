# MetaCompassv1.0
Last updated: July 12th, 2016



# REQUIREMENTS:
GNU C/C++; Perl 3; BLAST 2.2.31; Bowtie 2.2.9; BWA; kmer-mask
#todo:add software


# INSTALLATION:
    git clone https://github.com/marbl/MetaCompass.git
    cd MetaCompass
    wget https://gembox.cbcb.umd.edu/metacompass/refseq.tar.gz
    tar -xzvf refseq.tar.gz
    ./install.pl


# USAGE <under construction>

-- I have a set of metagenomic reads, and want to perform comparative assembly.

    snakemake <input FASTA/FASTQ reads> [options] 


-- I know the reference genomes, or I want to perform comparative assembly for a particular genome.

    snakemake -s 


-- You can try MetaCompass on a test data set.

For a big data set assuming reference genomes are known:

    snakemake -s

For a big data set assuming reference genomes are NOT known:

    snakemake -s

# Output:
      Assembled contigs:
                contigs_[iteration #].fasta
      
      Taxonomic profiling output from MetaPhyler:
                mc.classification; mc.[species|genus|family|order|class|phylum].taxprof
      
      Selected Reference genomes:
                mc.refseq.fasta
