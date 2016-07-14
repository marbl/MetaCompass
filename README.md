# MetaCompassv1.0
Last updated: July 13th, 2016



# REQUIREMENTS:
GNU C/C++; Perl 3; BLAST 2.4.0; Bowtie 2.2.9; BWA 0.7; kmer-mask
#todo:add software


# INSTALLATION:
    git clone https://github.com/marbl/MetaCompass.git
    cd MetaCompass
    ./install.sh


# USAGE
-- Go to the snakemake folder.

    cd snakemake
    
--Customize the configuration files (config.json and config2.json) as necessary. In the following example, we want to assembly 2 datasets, 'Sample1.fastq' and 'Sample2.fastq', and the data is located in the folder 'test' in the parent directory:

    {
        "reads": {
        "S1": ["Sample1", "Sample2"]
        },
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
