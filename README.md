# MetaCompass v1.0
Last updated: August 4th, 2017

# Required software:

* Python3 (>=) 3.1: https://www.python.org/download/releases/3.0/
* snakemake 3.7.1: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
* BLAST+ (>=) 2.4.0: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
* bowtie2  (>=) 2.2.9: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/ 
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

    python3 go_metacompass.py -P [read1.fq,read2.fq] -l [max read length]-o [output_folder] -m [min coverage] -t [ncpu]

-- I know the reference genomes, or I want to perform comparative assembly for a particular genome.

    python3 go_metacompass.py -r [references.fasta] -P [read1.fq,read2.fq] -o [output_folder] -m [min coverage] -t [ncpu]


# OUTPUT

--Output folder contains the following files:
    
    1) metacompass_output folder:
         Assembled contigs:
                metacompass_output/metacompass.final.ctg.fa
         Selected Reference genomes sequences:
                metacompass_output/metacompass.recruited.fa  
         Selected Reference genomes ids and taxonomy ids
                metacompass_output/metacompass.recruited.ids
    2) Assembly information
         metacompass.tsv file with the following information: 
            -contig ID
            -contig size
            -reference genome used (only reference-guided assembled contigs)
            -start and end position in the genome where contig originates from
            -taxonomic label of the genome
            -name of the genome  

                
# EXAMPLES

# Reference-guided assembly with known reference genomes (no reference selection).
-- Input data is available in the tutorial folder:

    Reference genome file:  Candidatus_Carsonella_ruddii_HT_Thao2000.fasta
    Metagenomic reads:      thao2000.1.fq
                            thao2000.2.fq	
-- Run:
   
     python3 go_metacompass.py --r tutorial/Candidatus_Carsonella_ruddii_HT_Thao2000.fasta -P tutorial/thao2000.1.fq,tutorial/thao2000.2.fq -o example1_output -m 3 -t 4

# Reference-guided assembly with reference selection.

-- Download and extract metagenomic sample:

    ftp://public-ftp.hmpdacc.org/Illumina/posterior_fornix/SRS044742.tar.bz2
   
    SRS044742/
        SRS044742.denovo_duplicates_marked.trimmed.1.fastq
        SRS044742.denovo_duplicates_marked.trimmed.2.fastq
        SRS044742.denovo_duplicates_marked.trimmed.singleton.fastq
-- Run:
   
     python3 go_metacompass.py -P SRS044742/SRS044742.denovo_duplicates_marked.trimmed.1.fastq,SRS044742/SRS044742.denovo_duplicates_marked.trimmed.2.fastq -U SRS044742/SRS044742.denovo_duplicates_marked.trimmed.singleton.fastq -o example2_output

  
Contact:
vcepeda{at}cs.umd.edu
