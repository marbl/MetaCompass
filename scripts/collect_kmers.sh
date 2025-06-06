#!/bin/bash

usage="$(basename "$0") [-f | --forward] path_to_forward_read \n\
                        [-r | --reverse] path_to_reverse_read \n\
                        [-u | --unpaired] path_to_unpaired_read \n\
                        [-refs | --references] path to reference file directory (Directory structure assumes NCBI dataset downloaded structure) \n\
                        [-ms | --kmer_size] kmer size \n\
                        [-o | --out]   path to output directory \n\
                        [-l | --log]   path to log file\n\
                        [-h | --help]  print out help message\n"

forward=""
reverse=""
unpaired=""
references=""
out=""

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

        -refs|--references)
            references="$2";
            shift;;
        
        -ms|--kmer_size)
            ms="$2";
            shift;;

        -o|--out)
            out="$2";
            shift;;

        -l|--log)
            log="$2";
            shift;;

        -h|--help)
            echo $usage;
            exit 0;;
    
    esac
    shift
done


if [[ $out = "" ]] || [[ ! -d "${out}" ]]; then
    echo "ERROR: Must provide valid output folder. You provided: ${out}" >> $log
    exit 1
fi


if [[ $references = "" ]] || [[ ! -d "${references}" ]]; then
    echo "ERROR: Must provide valid references folder. You provided: ${references}" >> $log
    exit 1
fi

check_valid_reads() {
    read_to_check=$1
    if [[ $read_to_check != "" ]] && [[ ! -f $read_to_check ]]; then
        echo "ERROR: Reads were provided but not valid: $read_to_check" >> $log
        exit 1
    fi
}

check_valid_reads "$forward"
check_valid_reads "$reverse"
check_valid_reads "$unpaired"

read_kmers=${out}"/read_kmers"
ref_kmers=${out}"/ref_kmers"

mkdir $read_kmers
mkdir $ref_kmers

# collect read kmers
if [[ -f "$forward" ]] && [[ -f "$reverse" ]] && [[ -f "$unpaired" ]]; then

    jellyfish count -F 3 -m $ms -s 100M -t 8 -o read_kmers.jf $forward $reverse $unpaired

elif [[ -f "$forward" ]] && [[ -f "$reverse" ]]; then

    jellyfish count -F 2 -m $ms -s 100M -t 8 -o read_kmers.jf $forward $reverse

elif [[ -f "$unpaired" ]]; then

    jellyfish count -m $ms -s 100M -t 8 -o read_kmers.jf $unpaired

else
    echo "ERROR: Input reads not found!" >> $log
    exit 1

fi

jellyfish dump read_kmers.jf > ${read_kmers}/read_kmers.fasta
rm read_kmers.jf

#collect references kmers
for reference in $(ls $references); do

    if [[ ! -d ${references}/${reference} ]]; then
        continue;
    fi

    fna=$(ls ${references}/${reference} | grep .fna)
    jellyfish count -m $ms -s 100M -t 8 -o ${reference}.jf ${references}/${reference}/${fna}
    jellyfish dump ${reference}.jf > ${ref_kmers}/${reference}.fasta
    rm ${reference}.jf
    
done

#check if no reference presents (output kmer folder for references is empty)
if [[ ! $(ls -A $ref_kmers) ]]; then
    echo "ERROR: Folder of references is empty or not in correct format." >> $log
    exit 1
fi
