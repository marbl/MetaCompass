# MetaCompass v1.0
Last updated: May 5th, 2017

# Required software:

* Python3 (>=) 3.1: https://www.python.org/download/releases/3.0/
* snakemake 3.7.1: https://bitbucket.org/snakemake/snakemake/src
* BLAST+ (>=) 2.4.0: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
* bowtie2  (>=) 2.2.4: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/ 
* kmer-mask (May 13th, 2015): https://sourceforge.net/p/kmer/code/HEAD/tree/trunk/
* samtools (>=) 1.x: http://samtools.sourceforge.net/ 
* MEGAHIT (>=) 1.0.6: https://github.com/voutcn/megahit
* java
* perl

# Memory and Disk Space Requirements.
You must have at least XXMB of memory and XXMB of hard disk space to perform a normal installation.

# INSTALLATION:

    git clone https://github.com/marbl/MetaCompass.git
    cd MetaCompass
    ./install.sh

# USAGE    

-- I have a set of metagenomic reads, and want to perform comparative assembly.

    python go_metacompass.py -P read1,read2 -o output -t $ncpu

-- I know the reference genomes, or I want to perform comparative assembly for a particular genome.

    python go_metacompass.py


-- You can try MetaCompass on a test data set.

For a mock data set assuming reference genomes are known:

    python go_metacompass.py

For a mock data set assuming reference genomes are NOT known:

    python go_metacompass.py


# Output:
      Assembled contigs:
                contigs_[iteration #].fasta
      
      Selected Reference genomes:
                mc.refseq.fasta
