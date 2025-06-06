#!/bin/bash

usage="$(basename "$0") [-f | --forward] path_to_forward_read \n\
                        [-r | --reverse] path_to_reverse_read \n\
                        [-u | --unpaired] path_to_unpaired_read \n\
                        [-m | --minlen] Min contig length for MEGAHIT \n\
                        [-t | --threads] number of threads \n\
                        [-o | --output] path_to_output_folder \n\
                        [-l | --log] path_to_log_file \n"

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

        -t|--threads)
            threads="$2";
            shift;;

        -m|--minlen)
            minlen="$2";
            shift;;
        
        -o|--output)
            output="$2";
            shift;;
        
        -l|--log)
            log="$2";
            shift;;
        
        *)
            echo $usage;
            exit 1
    esac
    shift
done

# Check input conditions
if [[ ! -f ${forward} && ! -f ${reverse} && ! -f ${unpaired} ]]; then
    echo "ERROR: Must provide either forward, reverse or unpaired reads. Got: ${forward}, ${reverse}, ${unpaired}."
    exit 1
fi

# MEGAHIT requires that the output directory does not exist
# if [[ ! -d "${output}" ]]; then
#     echo "ERROR: Must provide a valid output folder (-o flag). You provided: ${output}"
#     exit 1
# fi

if [[ $log = "" ]]; then
    echo "ERROR: Must provide a valid log. You provided: ${log}"
    exit 1
fi

echo "********** START RUNNING MEGAHIT **********" &>> $log # dreate log file / overwrite if it exists

cmd="time megahit -o $output --min-count 3 --min-contig-len $minlen -t $threads "

printf %s "-----Running denovo_assembly-----\n$(date)\n\n" >> $log

if [[ -f $forward && -f $reverse ]]; then
    cmd+="-1 $forward -2 $reverse "
fi

if [[ -f $unpaired ]]; then
    cmd+="-r $unpaired"
fi

# run MEGAHIT
{ $cmd 2> $log ; } 2>> $log

echo "" &>> $log
echo "********** FINISH RUNNING MEGAHIT **********" &>> $log