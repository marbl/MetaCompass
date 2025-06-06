#!/bin/bash

usage="$(basename "$0") [-f | --forward] path_to_forward_read \n\
                        [-r | --reverse] path_to_reverse_read \n\
                        [-u | --unpaired] path_to_unpaired_read \n\
                        [-i | --inputs] path_to_filtered_refs\n\
                        [-ref | --reference_db] path_to_reference_db \n\
                        [-fil | --filter_refs] filtering step was done or not (boolean) \n\
                                                + true: use filtered references from --inputs \n\
                                                + false: use references from database \n\
                        [-d | --deapth_of_coverage] deapth of coverage cutoff (int) \n\
                        [-b | --breadth_of_coverage] breadth of coverage (float) \n\
                        [-p | --percent_markers_covered] percent of markers covered (int) \n\
                        [-t | --threads] number of threads \n\
                        [-s | --scripts] scripts folder \n\
                        [-o | --output] path_to_output_folder \n\
                        [-l | --log] path_to_log_file \n\
                        [-h | --help] display help usage"

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

        -i|--inputs)
            inputs="$2";
            shift;;

        -ref|--reference_db)
            reference_db="$2";
            shift;;

        -fil|--filter_refs)
            filter_refs="$2";
            shift;;

        -d|--depth_of_coverage)
            deapth_of_coverage="$2";
            shift;;

        -b|--breadth_of_coverage)
            breadth_of_coverage="$2";
            shift;;
        
        -p|--percent_markers_covered)
            percent_markers_covered="$2";
            shift;;

        -t|--threads)
            threads="$2";
            shift;;

        -s|--scripts)
            scripts="$2";
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

if [[ ! -d "${output}" ]]; then
    echo "ERROR: Invalid output folder: ${output}"
    exit 1
fi

if [[ ! -d "${reference_db}/marker_index" ]]; then
    echo "ERROR: Invalid database"
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


if [[ ! -f "${scripts}/identify_candidate_genomes.py" ]] || \
   [[ ! -f "${scripts}/marker_coverage_per_candidate_genome.py" ]]; then
    
    echo "ERROR: Invalid scripts folder."
    echo "       Scripts folder must contain identify_candidate_genomes.py and marker_coverage_per_candidate_genome.py"
    exit 1
fi

align_reads(){

    gene=$1

    # 1. if filtering step was done but filtered ref is invalid
    if [[ $filter_refs = true ]] && [[ ! -f "$inputs/$gene.keep.fastq" ]]; then
        echo "ERROR: Invalid path to filtered references (flag: --inputs)"
        exit 1

    # 2. if filtering step was done and fitered ref is valid
    elif [[ $filter_refs = true ]] && [[ -f "$inputs/$gene.keep.fastq" ]]; then
      
        bowtie2 --threads $threads \
            -x ${reference_db}/marker_index/${gene}/${gene}_clustered \
            -U ${inputs}/$gene.keep.fastq | \
            samtools view -bS > ${output}/${gene}.match.bam

    #3. otherwise, use reference candidates from database
    else


        if [[ -f $forward ]] && [[ -f $reverse ]]; then
            
            minimap2 -ax sr -t $threads ${reference_db}/marker_index/${gene}/${gene}_clustered.fna ${output}/sorted_mapped_reads.fq | \
            samtools view -b > ${output}/${gene}.match.combined.bam
        fi

        if [[ -f $unpaired ]]; then
            bowtie2 --threads $threads \
            -x ${reference_db}/marker_index/${gene}/${gene}_clustered \
            -U $unpaired | \
            samtools view -bS > ${output}/${gene}.match.u.bam
        fi

        if [[ -f $forward ]] && [[ -f $reverse ]] && [[ -f $unpaired ]]; then
            #merge output files
            samtools merge ${output}/${gene}.match.bam \
                            ${output}/${gene}.match.1.bam \
                            ${output}/${gene}.match.2.bam \
                            ${output}/${gene}.match.u.bam
            rm ${output}/${gene}.match.1.bam
            rm ${output}/${gene}.match.2.bam
            rm ${output}/${gene}.match.u.bam
        
        elif [[ -f $forward ]] && [[ -f $reverse ]]; then
            #merge output files
            mv  ${output}/${gene}.match.combined.bam ${output}/${gene}.match.bam
        
        elif [[ -f $unpaired ]]; then
            mv ${output}/${gene}.match.u.bam ${output}/${gene}.match.bam #rename output file
        fi
    fi
}
SECONDS=0

combined_path=${reference_db}/marker_index/marker_clustered.fna
if [[ -f $forward ]] && [[ -f $reverse ]]; then
    minimap2 -ax sr -t $threads  $combined_path $forward $reverse | \
    samtools view -bS -F 4 - | \
    samtools bam2fq - | \
    seqkit sort - > ${output}/sorted_mapped_reads.fq
fi

for gene in $(ls ${reference_db}/marker_index); do
    SECONDS=0

    printf "%s\n" "----- Clustering for $gene -----" >> $log

    # align reads to marker gene reference
    align_reads $gene

    run_time=$(($SECONDS / 60))
    printf "%s\n" "Bowtie2 alignment took: $run_time min" >> $log

    # sort alignment file
    samtools sort ${output}/${gene}.match.bam -o ${output}/${gene}.match.sorted.bam
    rm ${output}/${gene}.match.bam

    # make temp file with marker gene lengths
    tail -n+2 ${reference_db}/marker_index/${gene}/${gene}_genome.tsv | cut -f 1,4 > ${output}/${gene}_genome.txt

    # calculate coverage
    # -max flag says that anything equal to or above the specified coverage is output into a single line
    bedtools genomecov -ibam ${output}/${gene}.match.sorted.bam \
                        -g ${output}/${gene}_genome.txt \
                        -max ${deapth_of_coverage} > ${output}/${gene}_genomeCov.txt

    run_time=$(($SECONDS / 60))
    printf "%s\n" "bedtools genomecov took: $run_time min" >> $log

    # Remove temp file with the marker gene lengths
    rm ${output}/${gene}_genome.txt

    # Extract markers with good breadth of coverage at specified depth
    # First column is the marker gene name, 
    # Second column is depth of coverage and 
    # fifth column is the percentage of the gene covered at the specified depth (breadth)
    echo -e "seq\tcoverage" > ${output}/${gene}_marker_cov.txt
    awk -v doc=${deapth_of_coverage} \
        -v boc=${breadth_of_coverage} \
        '{if($2 == doc && $5 >= boc){print $1"\t"$5}}' ${output}/${gene}_genomeCov.txt >> \
        ${output}/${gene}_marker_cov.txt

    run_time=$(($SECONDS / 60))
    printf "%s\n" "Awk took: $run_time min" >> $log

    # Markers are representative of a cluster of sequences
    # For each marker that is covered, pull out the corresponding sequences in the cluster
    python3 ${scripts}/marker_coverage_per_candidate_genome.py ${gene} \
    ${output}/${gene}_marker_cov.txt \
    ${output}/${gene}_marker_cov_per_genome.txt \
    ${reference_db}/marker_index

    run_time=$(($SECONDS / 60))
    printf "%s\n" "marker_coverage_per_candidate_genome.py took: $run_time min" >> $log

    # Remove unnecessary files
    rm ${output}/${gene}.match.sorted.bam
    
    # keep for now to add info to report
    mkdir ${output}/${gene}
    mv ${output}/${gene}_genomeCov.txt ${output}/${gene}
    mv ${output}/${gene}_marker_cov.txt ${output}/${gene}
    # rm ${output}/${gene}_genomeCov.txt
    # rm ${output}/${gene}_marker_cov.txt

    #rm $inputs/$gene.keep.*.fastq
done

# Concatenate coverage for all genes
cat ${output}/*_marker_cov_per_genome.txt >> ${output}/marker_cov_per_genome.txt
# rm ${output}/*_marker_cov_per_genome.txt 

# Output a table with the reference genomes and how many marker genes are covered
# Output a list of candidate reference genomes
{ time python3 ${scripts}/identify_candidate_genomes.py ${percent_markers_covered} ${output} ${reference_db} || echo "No suitable candidate genomes" 2>> $log ; } 2>> $log
rm ${output}/marker_cov_per_genome.txt
