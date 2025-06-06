#!/bin/bash

usage="$(basename "$0") [-f | --forward] path_to_forward_read \n\
                        [-r | --reverse] path_to_reverse_read \n\
                        [-u | --unpaired] path_to_unpaired_read \n\
                        [-ref | --reference_db] path_to_referene_db \n\
                        [-ms | --ms] kmer-mask kmer size \n\
                        [-c | --clean] kmer-mask clean percentage \n\
                        [-m | --match] kmer-mask match percentage \n\
                        [-rn | --readlen] kmer-mask kmer read length \n\
                        [-o | --output] path_to_output_folder \n\
                        [-l | --log] path_to_log_file \n\
                        [-h | --help] display help usage"

# set defaults
forward=""
reverse=""
unpaired=""

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

        -ref|--reference_db)
            reference_db=$2;
            shift;;

        -ms|--ms)
            ms="$2";
            shift;;
        
        -c|--clean)
            clean="$2";
            shift;;

        -m|--match)
            match="$2";
            shift;;

        -mask|--masking)
            masking="$2";
            shift;;
        
        -rn|--readlen)
            readlen="$2";
            shift;;

        -o|--output)
            output="$2";
            shift;;
        
        -l|--log)
            log="$2";
            shift;;
        
        -h|--help)
            echo $usage;
            exit 0;;

        *)
            echo $usage;
            exit 1;;
    esac
    shift
done

if [[ ! -d "${reference_db}/marker_index" ]]; then
    echo "ERROR: Invalid database"
    exit 1
fi

if [[ -f $forward && -f $reverse ]] || [[ -f $unpaired ]]; then 
    :
else 
    echo "ERROR: Either paired or unpaired reads must be valid"
    exit 1
fi

exit_on_error() {
    exit_code=$1
    error_message=$2
    if [[ $exit_code -ne 0 ]]; then
        echo $error_message
        exit $exit_code
    fi
}

check_valid_reads() {
    read_to_check=$1
    if [[ $read_to_check != "" ]] && [[ ! -f $read_to_check ]]; then
        echo "ERROR: Reads were provided but not valid: $read_to_check"
        exit 1
    fi
}

check_valid_reads "$forward"
check_valid_reads "$reverse"
check_valid_reads "$unpaired"

filter () {
    for gene in $(ls ${reference_db}/marker_index); do
        
        printf "%s\n" "----- Filtering for $gene -----" >> $log

        # run kmer-mask to filter references
        cmd="time kmer-mask \
            -mdb ${reference_db}/marker_index/${gene}/${gene}_clustered.ms28 \
            $1 \
            -ms $ms \
            -clean $clean \
            -match $match \
            -nomasking \
            -t 12 \
            -l $readlen \
            -o ${output}/${gene}"

        # TODO: Unsure if log needed. Make sure placement doesn't affect cluster_refs.sh.
        #echo "CMD: $cmd" >> $log
        { $cmd ; } 2>> $log
        
        #Exit if previous command failed
        exit_on_error $? "Last command executed: $cmd"

        # Keep only files containing masked sequences (found in the reference)
        cat ${output}/${gene}.match.*.fastq > ${output}/${gene}.keep.fastq
        rm ${output}/${gene}*clean*fastq
        rm ${output}/${gene}*murky*fastq
        rm ${output}/${gene}*mixed*fastq
        rm ${output}/${gene}*match*fastq

    done
    
}

if [[ $forward != "" ]] && [[ $reverse != "" ]]; then 
    filter "-1 $forward -2 $reverse "
fi

if [[ $unpaired != "" ]]; then
    filter "-1 $unpaired "
fi