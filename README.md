# MetaCompass v1.0
Last updated: May 8th, 2017

# Required software:

* Python3 (>=) 3.1: https://www.python.org/download/releases/3.0/
* snakemake 3.7.1: https://bitbucket.org/snakemake/snakemake/src
* BLAST+ (>=) 2.4.0: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
* bowtie2  (>=) 2.2.4: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/ 
* kmer-mask (May 13th, 2015): https://sourceforge.net/p/kmer/code/HEAD/tree/trunk/
* samtools (>=) 1.x: http://samtools.sourceforge.net/ 
* MEGAHIT (>=) 1.0.6: https://github.com/voutcn/megahit
* Java runtime (>=) 1.7 
* Perl5 (>=) 5.16

# Memory and Disk Space Requirements.
You must have at least 13GB of hard disk space to perform a normal installation.
You must have 8GB or more memory to allocate to the JVM (used by pilon).

# INSTALLATION:

    git clone https://github.com/marbl/MetaCompass.git
    cd MetaCompass
    ./install.sh

# USAGE    

-- I have a set of metagenomic reads, and want to perform reference-guided assembly.

    python go_metacompass.py -P [read1.fq,read2.fq] -o [output_folder] -t [ncpu]

-- I know the reference genomes, or I want to perform comparative assembly for a particular genome.

    python go_metacompass.py -r [references.fasta] -P [read1.fq,read2.fq] -o [output_folder] -t [ncpu]


# Output directory:
      Assembled contigs:
                metacompass_output/metacompass.final.ctf.fa
      
      Selected Reference genomes:
                metacompass_output/metacompass.recruited.fa
                
                
Contact:
vcepeda{at}cs.umd.edu
