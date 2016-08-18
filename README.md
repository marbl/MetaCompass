# MetaCompass v1.0
Last updated: August 9th, 2016

# Required software:

* BLAST+ (>=) 2.4.0: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
* bowtie2  (>=) 2.2.4: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/ 
* kmer-mask (May 13th, 2015): https://sourceforge.net/p/kmer/code/HEAD/tree/trunk/
* Python3 (>=) 3.1: https://www.python.org/download/releases/3.0/
* MEGAHIT (>=) 1.0.6: https://github.com/voutcn/megahit
* samtools (>=) 1.x: http://samtools.sourceforge.net/ 
* snakemake 3.7.1: https://bitbucket.org/snakemake/snakemake/src

# INSTALLATION:

    git clone https://github.com/marbl/MetaCompass.git
    cd MetaCompass
    ./install.sh

# USAGE
-- Go to the snakemake folder.

    cd snakemake
    
--Customize the configuration files (config.json and config2.json) as necessary. In the following example, we want to assemble one dataset, 'Sample1.fastq', with sample prefix of 'Sample1'

    {
        "reads": ["/path/to/reads/Sample1.fq"],
        "sample": Sample1,
        "reference":["../test/refseq.fna"],
        "prefix":"../test",
        "memory": 50,
        "nthreads": 12
    }

    

-- I have a set of metagenomic reads, and want to perform comparative assembly.

    snakemake

-- I know the reference genomes, or I want to perform comparative assembly for a particular genome.

    snakemake -s Snakefile_refgeno


-- You can try MetaCompass on a test data set.

For a mock data set assuming reference genomes are known:

    snakemake -s Snakefile_refgeno

For a mock data set assuming reference genomes are NOT known:

    snakemake

# Output:
      Assembled contigs:
                contigs_[iteration #].fasta
      
      Taxonomic profiling output from MetaPhyler:
                mc.classification; mc.[species|genus|family|order|class|phylum].taxprof
      
      Selected Reference genomes:
                mc.refseq.fasta
