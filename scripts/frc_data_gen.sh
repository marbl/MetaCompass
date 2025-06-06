#!/bin/bash

# Perform for ref-guided, denovo and merged contigs.
# 1. Prep reads into one file.
# 2. Create bowtie2 index out of contigs.
# 3. Align reads to contigs.
# 4. Run frc code.

usage="
Generate Feature Response Curve data for reference-guided, denovo and merged contigs.

bash $(basename "$0") -f <forward_reads> -r <reverse_reads> -u <unpaired_reads> \
                      -ref <ref-guided_contigs> -d <denovo_contigs> -m <merged_contigs> -o <out_dir>

where:
        [-f | --forward] path_to_forward_read.
        [-r | --reverse] path_to_reverse_read.
        [-u | --unpaired] path_to_unpaired_read.
        [-s | --split_reads] path to mapped and unmapped reads (Directory structure assumes output format from align_reads step of reference culling process).
        [-ref | --refguided_contigs] path to ref-guided contigs directory (Directory structure assumes output format from pilon step of ref-guided assembly process).
        [-d | --denovo_contigs] path to file containing denovo contigs.
        [-m | --merged_contigs] path to file containing merged contigs.
        [-frc | --frc] path to FRC executable.
        [-libavg | --libsize_avg] Average library size of reads.
        [-libstd | --libsize_sdv] Standard Deviation for library size of reads.
        [-oa | --orientation_a] Orientation of forward reads.
        [-ob | --orientation_b] Orientation of reverse reads.
        [-mpsize | --mp_size_z_cutoff] Mate pair size cutoff.
        [-mpedge | --mp_edge_z_cutoff] Mate edge size cutoff.
        [-doclogp | --doc_logp_cutoff] Depth of coverage cutoff (log(p)).
        [-o | --out]   path to output directory.
        [-l | --log]   path to log file.
        [-h | --help]  print usage."

forward=""
reverse=""
unpaired=""
split_reads=""
r_contigs=""
d_contigs=""
m_contigs=""
FRC=""
# Start default params for FRC code #
# data properties
libsize_avg="165"
libsize_sdv="92.8"
orientation_a="+"
orientation_b="+"
# statistical parameters
mp_size_z_cutoff="2.58"
mp_edge_z_cutoff="2.58"
doc_logp_cutoff="7"
# End default params for FRC code #
out=""
log=""
help=0

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

        -s|--split_reads)
            split_reads="$2";
            shift;;

        -ref|--refguided_contigs)
            r_contigs="$2";
            shift;;
        
        -d|--denovo_contigs)
            d_contigs="$2";
            shift;;

        -m|--merged_contigs)
            m_contigs="$2";
            shift;;

        -frc|--frc)
            FRC="$2";
            shift;;

        -libavg|--libsize_avg)
            libsize_avg="$2";
            shift;;

        -libstd|--libsize_sdv)
            libsize_sdv="$2";
            shift;;

        -oa|--orientation_a)
            orientation_a="$2";
            shift;;

        -ob|--orientation_b)
            orientation_b="$2";
            shift;;

        -mpsize|--mp_size_z_cutoff)
            mp_size_z_cutoff="$2";
            shift;;

        -mpedge|--mp_edge_z_cutoff)
            mp_edge_z_cutoff="$2";
            shift;;

        -doclogp|--doc_logp_cutoff)
            doc_logp_cutoff="$2";
            shift;;

        -o|--out)
            out="$2";
            shift;;

        -l|--log)
            log="$2";
            shift;;

        -h|--help)
            help=1;
            shift;;
    
    esac
    shift
done

if [ $help = 1 ]; then
    echo "$usage"
    exit 0
fi

#check reads
if [[ ("$forward" != '' && "$reverse" == '') || ("$forward" == '' && "$reverse" != '') ]]; then 
    echo "Error: Only have forward or reverse but need both."
    exit 1
elif [[ "$forward" == '' && "$reverse" == '' && "$unpaired" == '' ]]; then
    echo "Error: No reads provided."
    exit 1
fi

# only merged contigs required.
if [[ "$m_contigs" == '' ]]; then
    echo "Error: Need merged contigs file."
    exit 1
elif [[ ! -f "$m_contigs" ]]; then
    echo "Error: Merged contigs file path does not exist."
    exit 1
fi

if [[ "$out" == '' ]]; then
    echo "Error: Need output folder."
    exit 1
elif [[ ! -d "$out" ]]; then
    echo "Error: Output folder path does not exist."
    exit 1
fi
#check existence
if [[ "$forward" != '' && ! -f "$forward" ]]; then
    echo "Error: Forward reads path does not exist."
    exit 1
fi

if [[ "$reverse" != '' && ! -f "$reverse" ]]; then
    echo "Error: Reverse reads path does not exist."
    exit 1
fi

if [[ "$unpaired" != '' && ! -f "$unpaired" ]]; then
    echo "Error: Unpaired reads path does not exist."
    exit 1
fi

if [[ "$split_reads" != '' && ! -d "$split_reads" ]]; then
    echo "Error: Split reads reads path does not exist."
    exit 1
fi

if [[ "$FRC" != '' && ! -f "$FRC" ]]; then
    echo "Error: FRC executable path does not exist."
    exit 1
fi

get_data() {
    reads=$1
    contigs=$2
    out_path=$3
    name=$4

    sam=${out_path}"/"${name}".sam"
    bam=${out_path}"/"${name}".bam"
    bed_pre=${out_path}"/"${name}".pre.bed"
    bed_post=${out_path}"/"${name}".post.bed"
    fai=${contigs}".fai"
    out_file=${out_path}"/"${name}".tsv"

    # create index
    bowtie2_index=${out_path}"/idx/idx"
    bowtie2_output=${out_path}"/bowtie2_output.txt"

    mkdir $out_path
    mkdir ${out_path}"/idx"
    bowtie2-build $contigs $bowtie2_index > $bowtie2_output

    # create necessary files
    cat $reads | bowtie2 -U - -x $bowtie2_index -S $sam -p 8 >> $bowtie2_output
    samtools view -bS $sam > $bam
    samtools faidx $contigs
    bamToBed -i $bam > $bed_pre
    env LC_COLLATE=C sort -k 1,1 -k 4,4 -o $bed_post $bed_pre

    cmd="$FRC $bed_post $fai $out_file $libsize_avg $libsize_sdv $orientation_a $orientation_b $mp_size_z_cutoff $mp_edge_z_cutoff $doc_logp_cutoff"

    echo "Execute: $cmd" >> $log
    $cmd
}

if [[ -f $forward && -f $reverse && -f $unpaired ]]; then

    get_data "$forward $reverse $unpaired" "$m_contigs" "${out}/merged" "merged"

elif [[ -f $forward && -f $reverse ]]; then

    get_data "$forward $reverse" "$m_contigs" "${out}/merged" "merged"

elif [[ -f $unpaired ]]; then

    get_data "$unpaired" "$m_contigs" "${out}/merged" "merged"

fi

if [[ -f $forward && -f $reverse && -f $unpaired ]]; then

    unmapped_reads1=$(find $split_reads -name '*_unmapped.1*')
    unmapped_reads2=$(find $split_reads -name '*_unmapped.2*')
    unmapped_readsu=$(find $split_reads -name '*_unmapped.u*')

    get_data "$unmapped_reads1 $unmapped_reads2 $unmapped_readsu" "$d_contigs" "${out}/denovo" "denovo"

elif [[ -f $forward && -f $reverse ]]; then

    unmapped_reads1=$(find $split_reads -name '*_unmapped.1*')
    unmapped_reads2=$(find $split_reads -name '*_unmapped.2*')

    get_data "$unmapped_reads1 $unmapped_reads2" "$d_contigs" "${out}/denovo" "denovo"

elif [[ -f $unpaired ]]; then

    unmapped_readsu=$(find $split_reads -name '*_unmapped.u*')

    get_data "$unmapped_readsu" "$d_contigs" "${out}/denovo" "denovo"

fi

ref_out=$out"/ref_guided"
mkdir $ref_out

for ref in $(ls $split_reads); do

    out=${ref_out}"/"${ref}
    fasta="${r_contigs}/${ref}/contigs.pilon.fasta"

    if [[ ! -f $fasta ]]; then
        continue
    fi

    if [[ -f $forward && -f $reverse && -f $unpaired ]]; then

        mapped_reads1=$(find ${split_reads}/${ref} -name '*_mapped.1*')
        mapped_reads2=$(find ${split_reads}/${ref} -name '*_mapped.2*')
        mapped_readsu=$(find ${split_reads}/${ref} -name '*_mapped.u*')

        get_data "$mapped_reads1 $mapped_reads2 $mapped_readsu" "$fasta" "$out" "$ref"

    elif [[ -f $forward && -f $reverse ]]; then

        mapped_reads1=$(find ${split_reads}/${ref} -name '*_mapped.1*')
        mapped_reads2=$(find ${split_reads}/${ref} -name '*_mapped.2*')

        get_data "$mapped_reads1 $mapped_reads2" "$fasta" "$out" "$ref"

    elif [[ -f $unpaired ]]; then
    
        mapped_readsu=$(find ${split_reads}/${ref} -name '*_mapped.u*')

        get_data "$mapped_readsu" "$fasta" "$out" "$ref"

    fi
done