# MetaCompass v1.3
Last updated: January 16th, 2020


# Publication
Victoria Cepeda, Bo Liu, Mathieu Almeida, Christopher M. Hill, Sergey Koren, Todd J. Treangen, Mihai Pop.
bioRxiv 212506; doi: https://doi.org/10.1101/212506

# Required software:

* Python3 (>=) 3.1: https://www.python.org/download/releases/3.0/
* snakemake (>=) v3.7.1: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
* BLAST+ (>=) 2.4.0: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
* bowtie2  (>=) 2.2.9: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9
* kmer-mask (May 13th, 2015): https://sourceforge.net/projects/kmer/files/meryl-r2013.tar.xz
* mash (>=) 2.1: https://github.com/marbl/Mash/releases/tag/v2.1
* samtools (>=) 1.x: http://samtools.sourceforge.net/ 
* MEGAHIT (>=) 1.0.6: https://github.com/voutcn/megahit
* Java runtime (>=) 1.7 

# Memory and Disk Space Requirements.
You must have at least 90GB of hard disk space to perform a normal installation.
You must have 8GB or more memory to allocate to the JVM (used by pilon).

# INSTALLATION From Source (RECOMMENDED):
Get the Latest release from https://github.com/marbl/MetaCompass/releases:

    wget https://github.com/marbl/MetaCompass/archive/1.xx.tar.gz
    tar -xzvf 1.xx.tar.gz
    cd MetaCompass-1.xx
    ./install.sh

# INSTALLATION Using Git:

    git clone https://github.com/marbl/MetaCompass.git
    cd MetaCompass
    ./install.sh

# USAGE    

-- I have a set of metagenomic reads, and want to perform reference-guided assembly. 

    python3 go_metacompass.py -1 [read1.fq] -2 [read2.fq] -l [max read length]-o [output_folder] -m [min coverage] -t [ncpu] -y [memory]

-- I know the reference genomes, or I want to perform comparative assembly for a particular genome.

    python3 go_metacompass.py -r [references.fasta] -1 [read1.fq] -2 [read2.fq] -o [output_folder] -m [min coverage] -t [ncpu] -y [memory]


# OUTPUT

--Output folder contains the following files:
    
    1) metacompass_output folder:
         - Assembled contigs:
                metacompass_output/metacompass.final.ctg.fa
         - Selected Reference genomes sequences:
                metacompass_output/metacompass.recruited.fa  
         - Selected Reference genomes ids and taxonomy ids
                metacompass_output/metacompass.recruited.ids
    2) metacompass.tsv file with the following information: 
         - contig ID
         - contig size
         - reference genome used (only reference-guided assembled contigs)
         - name of the genome  
       

# EXAMPLES

# Reference-guided assembly with known reference genomes (no reference selection).
-- Input data is available in the tutorial folder:

    Reference genome file:  Candidatus_Carsonella_ruddii_HT_Thao2000.fasta
    Metagenomic reads:      thao2000.1.fq
                            thao2000.2.fq	
-- Run:
   
     python3 go_metacompass.py -r tutorial/Candidatus_Carsonella_ruddii_HT_Thao2000.fasta -1 tutorial/thao2000.1.fq -2 tutorial/thao2000.2.fq -l 150 -o example1_output -m 1 -t 4 -y 8

# Reference-guided assembly with reference selection.

-- Download and extract metagenomic sample:

    wget http://downloads.hmpdacc.org/dacc/hhs/genome/microbiome/wgs/analysis/hmwgsqc/v2/SRS044742.tar.bz2
    tar -xvf SRS044742.tar.bz2

-- The metagenomic sample contains:

    SRS044742/
        SRS044742.denovo_duplicates_marked.trimmed.1.fastq
        SRS044742.denovo_duplicates_marked.trimmed.2.fastq
        SRS044742.denovo_duplicates_marked.trimmed.singleton.fastq
-- Run:
   
     python3 go_metacompass.py -1 SRS044742/SRS044742.denovo_duplicates_marked.trimmed.1.fastq -2 SRS044742/SRS044742.denovo_duplicates_marked.trimmed.2.fastq -U SRS044742/SRS044742.denovo_duplicates_marked.trimmed.singleton.fastq -l 100 -o example2_output -y 8

  
Contact:
vcepeda@cs.umd.edu
