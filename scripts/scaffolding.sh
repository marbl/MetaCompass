#!/bin/bash

usage="$(basename "$0") \
                        [-f | --forward]            Forward forward.\n
                        [-r | --reverse]            Forward reverse.\n
                        [-u | --unpaired]           Forward unpaired.\n
                        [-c | --contigs]            File containings assembled contigs.\n
                        [-m | --metacarvel_path]    Path to MetaCarvel executable <run.py>.\n
                        [-o | --out]                Output directory.\n
                        [-l | --log]                Path to log file. Log file will be created, or overwritten if already existed \n
                                                        (default: <output_directory>/scffolding.log) \n
                        [-h | --help]               Print help message\n
                        
                        Takes in assembled contigs from the given reads and performs scaffolding\
                        by running MetaCarvel."


# set defaults
forward=""
reverse=""
unpaired=""
contigs=""
metacarvel_path=""

# Parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--forward)
            forward="$2";
            shift;;

        -r|--reverse)
            reverse="$2";
            shift;;

        -u|--unpaired)
            unpaired="$2";
            shift;;

        -c|--contigs)
            contigs="$2";
            shift;;

        -m|--metacarvel_path)
            metacarvel_path="$2";
            shift;;

        -o|--out)
            out="$2";
            shift;;
        
        -l|--log)
            log="$2";
            shift;;
        
        -h|--help)
            echo -e $usage;
            exit 0;;

        *)
            echo -e $usage;
            exit 1;;
    esac
    shift
done


# Check input conditions
if [[ ! -f ${forward} && ! -f ${reverse} && ! -f ${unpaired} ]]; then
    echo "ERROR: Must provide either forward, reverse or unpaired reads. Got: ${forward}, ${reverse}, ${unpaired}."
    exit 1
fi

if [[ ! -f ${contigs} ]]; then
    echo "ERROR: Must provide contigs. Got: ${contigs}"
    exit 1
fi

if [[ ! -f ${metacarvel_path} ]]; then
    echo "ERROR: Must provide path to MetaCarvel <run.py>. Got: ${metacarvel_path}"
    exit 1
fi

if [[ ! -d "${out}" ]]; then
    echo "ERROR: Must provide a valid output folder (-o flag). You provided: ${out}"
    exit 1
fi

if [[ $log = "" ]]; then
    log="${out}/valet.log"
fi

echo "********** START RUNNING METACARVEL **********" &> $log #Create log file / overwrite if it exists



# prepare date for metacarvel

idx=$out/idx
bowtie2-build $contigs $idx # build index
bam=""

if [[ -f $forward && -f $reverse ]]; then
    echo "Preparing data for $forward and $reverse." &> $log

    bowtie2 -x $idx -U $forward | samtools view -bS - | samtools sort - -o $out/alignment_1.bam
    bowtie2 -x $idx -U $reverse | samtools view -bS - | samtools sort - -o $out/alignment_2.bam

    bam+="$out/alignment_1.bam $out/alignment_2.bam "
fi

if [[ -f $unpaired ]]; then
    echo "Preparing data for $unpaired." &> $log

    bowtie2 -x $idx -U $unpaired | samtools view -bS - | samtools sort - -o $out/alignment_3.bam
    bam+="$out/alignment_3.bam"
fi

samtools merge $out/alignment_total.bam $bam
samtools sort -n $out/alignment_total.bam -o $out/alignment.bam

# run metacarvel

python3 $metacarvel_path -a $contigs -m $out/alignment.bam -d $out

echo "" &>> $log
echo "********** FINISH RUNNING METACARVEL **********" &>> $log