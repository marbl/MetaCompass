#!/bin/bash

#Print error message ($2) and exit with exit code ($1), if exit code is not zero
exit_on_error() {
    exit_code=$1
    error_message=$2
    if [[ $exit_code -ne 0 ]]; then
        echo "Exit code: $exit_code"
        echo $error_message
        exit $exit_code
    fi
}

run_command_and_exit_if_error() {
    cmd=$1
    $cmd #run cmd
    exit_code=$? #get exit code of $cmd

    #Parse exit code and error message to exit_on_error function
    exit_on_error $exit_code "Last command executed: $cmd"
}

# Run pilon for 1 reference
# Main output: <output>/<ref_name>/contigs.pilon
pilon_1_ref () {
    ref_name=$1
    ref_path=$2
    align_path=$3

    echo "Reference candidate: ${ref_path}" &>> $log
    echo "Sam file: ${align_path}" &>> $log

    # dont run if sam alignment file is empty
    if [[ $(wc -l < $align_path) == 0 ]]; then
        echo "Alignment file is empty. Not running pilon." &>> $log
        return
    fi

    #Create a subfolder inside output folder to store output for this reference candidate
    run_command_and_exit_if_error "mkdir ${output}/${ref_name}" &>> $log

    #Convert sam file to sorted bam file
    echo "" &>> $log
    echo "--- Converting sam file to sorted bam file ---" &>> $log
    run_command_and_exit_if_error "samtools view -bS ${align_path} \
                                    -o ${align_path}_out.sam" &>> $log

    run_command_and_exit_if_error "samtools sort -@ ${threads} ${align_path}_out.sam \
                                    -o ${align_path}_sorted.bam -O bam" &>> $log

    run_command_and_exit_if_error "samtools index ${align_path}_sorted.bam" &>> $log
    echo "--- Finish converting sam file to sorted bam file ---" &>> $log

    # Run pilon using sorted bam file
    echo "" &>> $log
    echo "--- Running pilon ---" &>> $log
    run_command_and_exit_if_error "java -Xmx${memory}G -jar $pilon \
                                    --flank ${flank} --threads ${threads} \
                                    --mindepth ${mindepth} --genome ${ref_path} \
                                    --frags ${align_path}_sorted.bam \
                                    --output contigs.pilon --outdir ${output}/${ref_name} \
                                    --fix bases,amb --tracks --changes" &>> $log
    
}

usage="$(basename "$0") [-al | --alignment] directory containing alignments for the reference candidates.\n\
                                            Each sub-directory corresponds to 1 reference candidate and has 1 alignment.sam file \n\
                                            Do not include / at the end of the directory\n\
                        [-ref | --reference] input reference folder, has the same format as the marker_index folder of NCBI database.\n\
                                            Do not include / at the end of the directory\n\
                        [-f | --flank] input flank arg to pilon (default=5)\n\
                        [-t | --threads] input threads arg to pilon (default=8)\n\
                        [-md | --mindepth] input mindepth arg to pilon (default=3) \n\
                        [-m | --memory] max memory (in GB) assigned for pilon call (default=20)\n\
                        [-pl | --pilon-path] path to pilon jar file\n\
                        [-o | --out] output directory for pilon.\n\
                                    Do not include / at the end of the directory \n\
                        [-l | --log] path to log file. Log file will be created, or overwritten if already existed \n\
                                    (default: <output_directory>/pilon.log) \n\ 
                        [-h | --help] print help message\n"


# set defaults
flank="5"
threads="12"
mindepth="3"
reference=""
output=""
alignment=""
memory="20"
pilon=""
log=""

# Parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -al|--alignment)
            alignment="$2";
            shift;;
        
        -f|--flank)
            flank="$2";
            shift;;
        
        -t|--threads)
            threads="$2";
            shift;;

        -md|--mindepth)
            mindepth="$2";
            shift;;
        
        -ref|--reference)
            reference=$2;
            shift;;
        
        -m|--memory)
            memory=$2;
            shift;;

        -pl|--pilon-path)
            pilon=$2;
            shift;;

        -o|--out)
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


# Check input conditions
if [[ ! -d "${output}" ]]; then
    echo "ERROR: Must provide a valid output folder (-o flag). You provided: ${output}"
    exit 1
fi

if [[ $log = "" ]]; then
    log="${output}/pilon.log"
fi

if [[ ! -f "${pilon}"  ]]; then
    echo "ERROR: Pilon path (-pl flag) must be provided and be valid. You provided: ${pilon}"
    exit 1
fi

if [[ ! -d "${alignment}"  ]]; then
    echo "ERROR: Directory containing alignments (-al flag) must be valid. You provided: ${alignment}"
    exit 1
fi

if [[ ! "$(ls ${alignment})" ]]; then
    echo "ERROR: Invalid or empty alignment folder (-al flag). You provided: ${alignment}"
    exit 1
fi

if [[ ! -d "${reference}" ]]; then
    echo "ERROR: Directory containing references (-ref flag) must be valid. You provided: ${reference}"
    exit 1
fi

echo "********** START RUNNING PILON FOR ALL REFERENCE CANDIDATES **********" &> $log #Create log file / overwrite if it exists

# For each reference candidate in the alignment directory (-al),
# find the corresponding reference fasta file in -ref folder, then run pilon
for gene in $(ls ${alignment}); do
    echo "" &>> $log
    echo "##### REFERENCE: ${gene} #####" &>> $log
    
    # Check if alignment.sam exists for this gene
    sam_file_paired="${alignment}/${gene}/paired/data/alignment.sam"
    sam_file_interleaved="${alignment}/${gene}/interleaved/data/alignment.sam"
    if [[ ! -f "${sam_file_paired}" && ! -f "${sam_file_interleaved}" ]]; then
        echo "ERROR: Cannot find alignment.sam in ${sam_file_paired} or ${sam_file_interleaved}" &>> $log
        continue
    fi
    
    # Check if this same gene appears in -ref folder
    if [[ ! -d "${reference}/${gene}" ]]; then
        echo "ERROR: Found ${gene} in ${alignment} but could not find ${gene} in ${reference}" &>> $log
        echo "      ${reference}/${gene} is not a valid folder" &>> $log
        exit 1
    fi

    # Get the reference fasta file from -ref folder.
    # Assume there is only 1 fasta file with pattern ${gene}_*.fna in ${reference}/${gene}
    ref_file=""
    for file in $(find "${reference}/${gene}" -type f -name "${gene}_*.fna"); do
        ref_file="${file}"
    done

    if [[ $ref_file = "" ]]; then
        # Grab any fna
        # TODO: Need to account for multiple fna types or look for any fna file.
        ref_file=$(ls ${reference}/${gene}/*.fna)

        if [[ $ref_file = "" ]]; then
            echo "ERROR: Reference fasta file with pattern ${gene}_*.fna cannot be found in ${reference}/${gene}" &>> $log
            exit 1
        fi
    fi

    #Run pilon for the current reference, $1: reference name, $2: path to reference, $3: path to alignment.sam file
    if [[ -f "${sam_file_paired}" ]]; then
        pilon_1_ref "${gene}" "${ref_file}" "${sam_file_paired}" &>> $log
    fi
    if [[ -f "${sam_file_interleaved}" ]]; then
        pilon_1_ref "${gene}" "${ref_file}" "${sam_file_interleaved}" &>> $log
    fi

    echo "##### FINISH RUNNING PILON FOR REFERENCE: ${gene} #####" &>> $log
    
done

echo "" &>> $log
echo "********** FINISH RUNNING PILON FOR ALL REFERENCE CANDIDATES **********" &>> $log