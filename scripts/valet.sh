#!/bin/bash

usage="$(basename "$0") \
                        [-r | --reads]              Directory containing folder names of references containing reads aligned to that particular reference.\n
                        [-orient | --orientation]   Orientation of reads.\n
                        [-refs | --references]      Directory to references.\n
                        [-valet | --valet]          Path to valet executable.\n
                        [-run | --run_valet]        Whether or not to run VALET. \n
                        [-vis | --visualize]        Path to visualize script for valet. Set only if visualizations desired.\n
                        [-o | --out]                Output directory for pilon.\n
                        [-l | --log]                Path to log file. Log file will be created, or overwritten if already existed \n
                                                        (default: <output_directory>/valet.log) \n
                        [-h | --help]               Print help message\n
                        
                        Takes in a directory containing directories with references and read files aligned to them.
                        Produces fai files for the referenes and then runs valet."


# set defaults
reads=""
orientation=""
references=""
valet=""
VISUALIZE=""


# Parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -r|--reads)
            reads="$2";
            shift;;
            
        -orient|--orientation)
            orientation="$2";
            shift;;

        -refs|--references)
            references="$2";
            shift;;

        -val|--valet)
            valet="$2";
            shift;;

        -run|--run_valet)
            run_valet="$2";
            shift;;

        -vis|--visualize)
            VISUALIZE="$2";
            shift;;

        -o|--out)
            output="$2";
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
if [[ ! -d "${output}" ]]; then
    echo "ERROR: Must provide a valid output folder (-o flag). You provided: ${output}"
    exit 1
fi

if [[ $log = "" ]]; then
    log="${output}/valet.log"
fi

if [[ "${valet}" = "" ]]; then
    echo "ERROR: Must provide path to velet executable. You provided: ${valet}"
    exit 1
fi

if [[ ! -d ${reads} ]]; then
    echo "ERROR: Must provide path to reads. You provided: ${reads}"
    exit 1
fi

if [[ ${orientation} = "" ]]; then
    echo "Did not provide orientation of reads. Defaulting to --fr"
    orientation="--fr"
fi

if [[ ! -d ${references} ]]; then
    echo "ERROR: Must provide path to references. You provided: ${references}"
    exit 1
fi

if [[ "${VISUALIZE}" != "" ]] && [[ ! -f "${VISUALIZE}" ]]; then
    echo "ERROR: Must provide correct script for visualization. You provided: ${VISUALIZE}"
    exit 1
fi

echo "********** START RUNNING VALETS FOR ALL REFERENCES **********" &> $log #Create log file / overwrite if it exists

for reference in $(ls $reads); do

    # get reads for valet
    ref_fasta=$(ls ${references}/${reference}/*.fna) # TODO could be .fq, .fasta, .fa, .ffa, etc..., but fai file will be in the same directory
    r1="${reads}/${reference}/reads_mapped.1.fq"
    r2="${reads}/${reference}/reads_mapped.2.fq"
    ru="${reads}/${reference}/reads_mapped.u.fq"

    out=${output}/${reference}
    cmd=""

    # create fai
    samtools faidx $ref_fasta
    fai=$(ls ${references}/${reference}/*.fai) # fai created where fasta file is located
    
    # paired reads and visualize if necesssary
    if [[ -f $r1 && $(wc -l < $r1) > 0 ]] && [[ -f $r2 && $(wc -l < $r2) > 0 ]]; then
        
        cmd="${valet} --fasta ${ref_fasta} --fai ${fai} --paired ${r1} ${r2} ${orientation} --output ${out}/paired"

        if [[ ${VISUALIZE} != "" ]]; then
            vis="${VISUALIZE} ${out}/paired ${fai} ${r1} ${r2}"
        fi

        # create output dir
        mkdir $out

        if [[ $run_valet == "true" ]]; then

            echo "Running: $cmd" >> $log
            $cmd

            if [[ ${VISUALIZE} != "" ]]; then
                echo "Running: $vis" >> $log
                $vis
            fi

        elif [[ $run_valet == "false" ]]; then

            out_dir=${out}/paired/data
            mkdir -p $out_dir
            index=${out_dir}/index
            mkdir $index

            # align mapped reads to reference and create sam file
            bowtie2-build $ref_fasta $index/index
            bowtie2 -x ${index}/index -1 ${r1} -2 ${r2} -S ${out_dir}/alignment.sam

        fi

    fi

    # interleaved reads and visualize if necessary
    if [[ -f $ru && $(wc -l < $ru) > 0 ]]; then

        cmd="${valet} --fasta ${ref_fasta} --fai ${fai} --interleaved ${ru} --output ${out}/interleaved"

        if [[ ${VISUALIZE} != "" ]]; then
            vis="${VISUALIZE} ${out}/interleaved ${fai} ${ru}"
        fi

        # create output dir
        mkdir $out

        if [[ $run_valet == "true" ]]; then

            echo "Running: $cmd" >> $log
            $cmd
            
            if [[ ${VISUALIZE} != "" ]]; then
                echo "Running: $vis" >> $log
                $vis
            fi

        elif [[ $run_valet == "false" ]]; then

            out_dir=${out}/interleaved/data
            mkdir -p $out_dir
            index=${out_dir}/index
            mkdir $index

            # align mapped reads to reference and create sam file
            bowtie2-build $ref_fasta $index/index
            bowtie2 -x ${index}/index -U ${ru} -S ${out_dir}/alignment.sam

        fi
    fi

done

echo "" &>> $log
echo "********** FINISH RUNNING VALET FOR ALL REFERENCES **********" &>> $log