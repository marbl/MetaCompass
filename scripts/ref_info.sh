#!/bin/bash

usage="Used to calculate number of sequences mapped to a reference, as well as its breadth and depth of coverage. \n\n\
        
        $(basename "$0")    [-f | --forward] Path to forward read. \n\
                            [-r | --reverse] Path to reverse read. \n\
                            [-u | --unpaired] Path unpaired reads. \n\
                            [-refs | --references] Path to folder containing references. \n\
                            [-d | --depth] Depth of coverage for depth and breadth. \n\
                            [-ru | --run] Skip calculations and fill with NA's. \n\
                            [-o | --output] Path to output folder. \n\
                            [-l | --log] Path to log file. \n\
                            [-h | --help] Display help usage."

forward=""
reverse=""
unpaired=""
references=""

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
        
        -d|--depth)
            depth="$2";
            shift;;

        -ru|--run)
            run="$2";
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
            echo "You did something wrong!!!"
            echo $usage;
            exit 1;;
    esac
    shift
done

if [[ ! -d ${output} ]]; then
    echo "ERROR: Invalid output folder: ${output}"
    exit 1
fi

if [[ ! -d "${references}" ]]; then
    echo "ERROR: Invalid references folder: ${references}"
    exit 1
fi

if [[ -f $forward && -f $reverse ]] || [[ -f $unpaired ]]; then 
    :
else 
    echo "ERROR: Either paired or unpaired reads must be valid"
    exit 1
fi

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

out_file="${output}/ref_metrics.tsv"
printf "Reference ID\tNum Reads Mapped\tDepth of Coverage\tBreadth of Coverage\n" >> $out_file

if [[ $run == "true" ]]; then
    for reference in $(ls $references); do

        if [[ ! -d ${references}/${reference} ]]; then
            continue
        fi

        # build index
        idx="${references}/${reference}/idx"
        bam="${references}/${reference}/mapping_result_sorted.bam"
        fna=$(ls ${references}/${reference} | grep ".fna$")

        if [[ ! -f "${references}/${reference}/${fna}" ]]; then
            printf "Could not get .fna file for $reference\n"
            printf "Got: $fna\n"
            continue
        fi
        
        bowtie2-build ${references}/${reference}/${fna} $idx

        # map reads and sort bam file
        if [[ -f $forward ]] && [[ -f $reverse ]]; then
            bowtie2 -x $idx --no-unal -1 $forward -2 $reverse -S - -p 12 | \
                samtools view -bS - | \
                samtools sort -m 5G >> $bam
        elif [[ -f $unpaired ]]; then
            bowtie2 -x $idx --no-unal -U $unpaired -S - -p 12 | \
                samtools view -bS - | \
                samtools sort -m 5G >> $bam
        fi

        samtools index $bam

        # calculations
        num_reads_covered=$(samtools view -F 0x4 $bam | cut -f 1 | sort | uniq | wc -l)

        num_bases_covered_at_d=$(samtools mpileup $bam | awk -v X="${depth}" '$4>=X' | wc -l)
        len_reference=$(bowtie2-inspect -s $idx | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}')

        breadth=$(python -c "print('%.2f' % ($num_bases_covered_at_d / $len_reference * 100))") # use python to do float division
        depth_sum=$(samtools depth -a $bam  |  awk '{sum+=$3} END { print sum }') # calc for covered and non-covered regions. remove "-a" for just covered regions
        depth=$(python -c "print('%.2f' % ($depth_sum / $len_reference))") # use python to do float division
        # write to file
        printf "%s\t%d\t%f\t%f\n" "$reference" "$num_reads_covered" "$depth" "$breadth" >> $out_file

    done
else

    for reference in $(ls $references); do

        if [[ ! -d ${references}/${reference} ]]; then
            continue
        fi

        num_reads_covered='NA'
        depth='NA'
        breadth='NA'

        # write to file
        printf "%s\t%s\t%s\t%s\n" "$reference" "$num_reads_covered" "$depth" "$breadth" >> $out_file
    done

fi