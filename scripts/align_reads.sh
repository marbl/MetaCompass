#!/bin/bash

N=`date +%s%N`
export PS4='+[$(((`date +%s%N`-$N)/1000000))ms][${BASH_SOURCE}:${LINENO}]: ${FUNCNAME[0]:+${FUNCNAME[0]}(): }';
exec 2>runtime.log; set -x;


usage="
Align reads to a list of references sequentially by aligning to first reference,
then taking unaligned reads and aligning to next reference until no reads or references left.

bash $(basename "$0") -f <forward_reads> -r <reverse_reads> -u <unpaired_reads> -ref <reference_dir> -min-list <min_list_file> -o <out_dir>

where:
        [-f | --forward] path_to_forward_read
        [-r | --reverse] path_to_reverse_read
        [-u | --unpaired] path_to_unpaired_read
        [-refs | --references] path to reference file directory (Directory structure assumes NCBI dataset downloaded structure)
        [-min-list | --minimum-list] path to file containing minimum list of references to use, in order (descending)
        [-o | --out]   path to output directory
        [-l | --log]   path to log file
        [-p | --pilon-path]   path to pilon jar file (metacompass/bin/pilon-1.23.jar)
        [-t | --threads] input threads arg to bowtie2 (default=8)\n\
        [-s | --stat-script] path to stat-script
	[-h | --help]  print usage"

forward=""
reverse=""
unpaired=""
references=""
min_list=""
stat_script=""
pilon_path=""
out=""
help=0
threads=8  
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

        -refs|--references)
            references="$2";
            shift;;
        
        -min-list|--minimum-list)
            min_list="$2";
            shift;;

        -o|--out)
            out="$2";
            shift;;

	-s|--stat-script)
            stat_script="$2";
            shift;;	


        -l|--log)
            log="$2";
            shift;;

        -p|--pilon-path)
            pilon_path="$2";
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
#check for other required fields
if [[ "$references" == '' ]]; then
    echo "Error: Need references."
    exit 1
elif [[ ! -d "$references" ]]; then
    echo "Error: References path does not exist."
    exit 1
fi

if [[ "$min_list" == '' ]]; then
    echo "Error: Need min-list file."
    exit 1
elif [[ ! -f "$min_list" ]]; then
    echo "Error: Min-list file path does not exist."
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

elif [[ "$forward" != '' && -f "$forward" ]]; then

    tmp_num=$(wc -l < $forward)
    num_reads=$((tmp_num/4 | bc))
fi

if [[ "$reverse" != '' && ! -f "$reverse" ]]; then
    echo "Error: Reverse reads path does not exist."
    exit 1
fi

if [[ "$unpaired" != '' && ! -f "$unpaired" ]]; then
    echo "Error: Unpaired reads path does not exist."
    exit 1

elif [[ "$unpaired" != '' && -f "$unpaired" ]]; then

    tmp_num=$(wc -l < $unpaired)
    num_reads=$((tmp_num/4 | bc))
fi

if [[ "$log" != '' && ! -f "$log" ]]; then
    echo "Error: Log file path does not exist."
    exit 1
fi

run=0 # track if there are reads left to align
num_mapped_reads=0
check=0 # Track if input reads split. Need to remove subsequent unmapped reads to prevent excessive space usage.

MIN_COVERAGE_DEPTH=1
# read references
while IFS= read -r reference || [[ -n "$reference" ]]
do
    # reference=$(echo $reference | cut -d '.' -f1)".1" 

    # create bowtie2 index and output dir for reference
    out_dir=${out}/$reference
    mkdir ${out_dir}
    index=${out_dir}/index
    mkdir $index
    # TODO: Handle different file extensions
    ref_file=$(ls ${references}/$reference | grep ".fna$" | grep -v "cds_from_genomic.fna" )
    echo $ref_file 
    if [[ ! -f ${references}/$reference/$ref_file ]]; then

        printf "%s\n" "Metamorphic relation potentially violated: Minimum reference set is not subset of all references selected!!!" >> $log
        printf "%s\n" "Reference .fna file does not exist for $reference!!!" >> $log
        printf "%s\n" "Reference file: $ref_file." >> $log
        printf "%s\n" "Please check parameters and contact support if unresolved!!!" >> $log
        exit 1
    fi

    bowtie2-build  ${references}/$reference/$ref_file $index/index
    nb_ref=$(bowtie2-inspect -s  $index/index | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}')
    echo "check1"
    wc -l $forward
    wc -l $reverse
    echo $unpaired
   
    # Alignment
    if [[ -s $forward &&  -s $reverse && $(wc -l < $forward) > 0 && $(wc -l < $reverse) > 0 ]]; then
        echo "check2"
	if [[ -s $unpaired && $(wc -l < $unpaired) > 0 ]]; then
           bowtie2 -x $index/index -1 $forward -2 $reverse -U $unpaired -S ${out_dir}/alignment.sam -p ${threads}
        else

           bowtie2 -x $index/index -1 $forward -2 $reverse -S ${out_dir}/alignment.sam -p ${threads}
        fi

	
        samtools view --threads ${threads} -bS -h ${out_dir}/alignment.sam > ${out_dir}/alignment_temp.bam  
	samtools view --threads ${threads} -h -F 4  ${out_dir}/alignment_temp.bam >  ${out_dir}/alignment_temp_filt.bam 
	rm ${out_dir}/alignment_temp.bam
	samtools sort -@ ${threads} ${out_dir}/alignment_temp_filt.bam  -o ${out_dir}/alignment_sorted.bam
	rm  ${out_dir}/alignment_temp_filt.bam
	samtools index ${out_dir}/alignment_sorted.bam
	samtools depth ${out_dir}/alignment_sorted.bam -J  > ${out_dir}/depth.txt
        awk '{print $1"\t"($2)-1"\t"$2"\t"$3}'  ${out_dir}/depth.txt> ${out_dir}/non_zero_region.bed
        rm ${out_dir}/depth.txt
        bedtools merge -i  ${out_dir}/non_zero_region.bed |awk '{print $1":"($2)+1"-"$3}' > ${out_dir}/cord.txt
	samtools faidx ${references}/$reference/$ref_file -r  ${out_dir}/cord.txt > ${out_dir}/non_zero_region_contig.fasta
	sed -i 's/:/_/g ; s/-/_/g ;  s/\./_/g' ${out_dir}/non_zero_region_contig.fasta
        seqkit seq -m 500 ${out_dir}/non_zero_region_contig.fasta > ${out_dir}/non_zero_region_contig_filt.fasta

    if [[ $(wc -l < ${out_dir}/non_zero_region_contig_filt.fasta) -eq 0 ]]; then 
	    continue
	fi 

	rm ${out_dir}/non_zero_region_contig.fasta 
        mv ${out_dir}/non_zero_region_contig_filt.fasta ${out_dir}/non_zero_region_contig.fasta 
        samtools view --threads ${threads} -h -F 8 -f 0x40 ${out_dir}/alignment_sorted.bam | samtools  bam2fq > ${out_dir}/reads_mapped.1.fq
        samtools view --threads ${threads} -h -F 8 -f 0x80 ${out_dir}/alignment_sorted.bam | samtools  bam2fq >  ${out_dir}/reads_mapped.2.fq  
        samtools view --threads ${threads} -h -f 8  ${out_dir}/alignment_sorted.bam | samtools bam2fq > ${out_dir}/reads_mapped.u.fq
        samtools view --threads ${threads} -h -f 4   -bS ${out_dir}/alignment.sam >  ${out_dir}/alignment_unmapped.bam 
        
        #alignment to contigs 
        bowtie2-build ${out_dir}/non_zero_region_contig.fasta ${out_dir}/non_zero_region_contig
         if [[ $(wc -l < ${out_dir}/reads_mapped.u.fq) > 0 ]]; then
            bowtie2 -x ${out_dir}/non_zero_region_contig -1 ${out_dir}/reads_mapped.1.fq -2 ${out_dir}/reads_mapped.2.fq  -U ${out_dir}/reads_mapped.u.fq -S ${out_dir}/non_zero_contig_alignment.sam -p ${threads}
        else

            bowtie2 -x ${out_dir}/non_zero_region_contig -1 ${out_dir}/reads_mapped.1.fq -2 ${out_dir}/reads_mapped.2.fq  -S ${out_dir}/non_zero_contig_alignment.sam -p ${threads}
        fi

       
        samtools view --threads ${threads}  -h -bS ${out_dir}/non_zero_contig_alignment.sam >  ${out_dir}/non_zero_contig_alignment_temp.sam 
	samtools view  --threads ${threads}  -h -F 4 ${out_dir}/non_zero_contig_alignment_temp.sam > ${out_dir}/non_zero_contig_alignment_filt.sam
	samtools sort  -@ ${threads} ${out_dir}/non_zero_contig_alignment_filt.sam  -o ${out_dir}/non_zero_contig_alignment_sorted.bam 
	samtools index ${out_dir}/non_zero_contig_alignment_sorted.bam
        samtools view  --threads ${threads} -h -f 4  -bS ${out_dir}/non_zero_contig_alignment_temp.sam > ${out_dir}/non_zero_contig_alignment_unmapped.bam
        samtools merge --threads ${threads}  -o ${out_dir}/ummapped.sam ${out_dir}/non_zero_contig_alignment_unmapped.bam  ${out_dir}/alignment_unmapped.bam 
        
        samtools view  --threads ${threads} -h -f 8 -f 0x40 ${out_dir}/ummapped.sam | samtools  bam2fq  > ${out_dir}/reads_unmapped.1.fq
        samtools view  --threads ${threads} -h -f 8 -f 0x80 ${out_dir}/ummapped.sam | samtools  bam2fq  > ${out_dir}/reads_unmapped.2.fq
        samtools view  --threads ${threads} -h -F 8  ${out_dir}/ummapped.sam | samtools bam2fq  > ${out_dir}/reads_unmapped.u.fq
        rm ${out_dir}/reads_mapped.1.fq
        rm ${out_dir}/reads_mapped.2.fq
        rm ${out_dir}/reads_mapped.u.fq
        samtools bam2fq ${out_dir}/non_zero_contig_alignment_filt.sam > ${out_dir}/reads_mapped.fq
        #rm ${out_dir}/non_zero_contig_alignment_filt.sam
        samtools depth ${out_dir}/non_zero_contig_alignment_sorted.bam -J |  awk '{arr[$1]+=$3;} END {for (i in arr) print i, arr[i]}' |  sed 's/_/\t/g '| awk '{if($5>0 ){print $1"_"$2"_"$3"_"$4}}'> ${out_dir}/seq_list.txt
        seqkit grep -f ${out_dir}/seq_list.txt ${out_dir}/non_zero_region_contig.fasta >${out_dir}/non_zero_region_contig_cleaned.fasta
        mkdir $out_dir/pilon
        #pilon 
        java -Xmx20G -jar ${pilon_path} --flank 5 --threads ${threads} --mindepth 3 --genome ${out_dir}/non_zero_region_contig_cleaned.fasta  --frags ${out_dir}/non_zero_contig_alignment_sorted.bam  --output contigs.pilon --outdir $out_dir/pilon --fix bases,amb --tracks --changes 
        if [[ $check == 1 ]]; then
        
            rm $forward
            rm $reverse
            rm $unpaired
        else
            check=1
        fi
 
	  
        forward=${out_dir}/reads_unmapped.1.fq 
        reverse=${out_dir}/reads_unmapped.2.fq
        unpaired=${out_dir}/reads_unmapped.u.fq
        
    
	sort -t'_' -k1,1 -k2,2 -k3,3n ${out_dir}/seq_list.txt > ${out_dir}/seq_list_sorted.txt      	
        python ${stat_script}/AGP_generator.py -o ${out_dir} -c  ${out_dir}/seq_list_sorted.txt
	rm ${out_dir}/seq_list_sorted.txt
	echo Average Depth of Coverage, >>  ${out_dir}/stats
	samtools depth  ${out_dir}/non_zero_contig_alignment_sorted.bam  -J  |  awk '{arr[$1]+=$3;} END {for (i in arr) print i, arr[i]}' |  sed 's/_/\t/g '| awk '{s=s+$4-$3+1;g=g+$5} END {print g/s}' >>${out_dir}/stats
	echo Average Depth of Coverage per contig, >>  ${out_dir}/stats
	bedtools genomecov -d -ibam  ${out_dir}/non_zero_contig_alignment_sorted.bam |  awk '{arr[$1]+=$3;} END {for (i in arr) print i, arr[i]}' |  sed 's/_/\t/g '| awk '{if($5!=0){print $1"_"$2"_"$3"_"$4"\t" $5/($4-$3+1)}}' >>${out_dir}/stats
	n50 --format tsv ${out_dir}/non_zero_region_contig_cleaned.fasta | sed -n '2 p' |awk '{print "#seqs:\n" $2 "\nTotal length:\n" $3 "\nMax:\n" $6 }' >>${out_dir}/stats
    seqkit fx2tab --length --name ${out_dir}/non_zero_region_contig_cleaned.fasta >  ${out_dir}/length_assembly.txt
    seqkit fx2tab --length --name ${references}/$reference/$ref_file >  ${out_dir}/length_ref.txt
	python ${stat_script}/n50.py -i  ${out_dir}/length_assembly.txt -r ${out_dir}/length_ref.txt >>${out_dir}/stats
    rm  ${out_dir}/length_ref.txt
    rm  ${out_dir}/length_assembly.txt
	tmp_num=$(wc -l < ${out_dir}/reads_mapped.1.fq)
        tmp_num=$((tmp_num/4 | bc))
        num_mapped_reads=$((num_mapped_reads + tmp_num))
        run=1
	rm ${out_dir}/non_zero_contig_alignment_temp.sam
	rm ${out_dir}/non_zero_contig_alignment_unmapped.bam
    rm ${out_dir}/alignment_unmapped.bam 
	rm ${out_dir}/non_zero_contig_alignment_filt.sam
    rm ${out_dir}/ummapped.sam
    rm ${out_dir}/seq_list.txt
    rm ${out_dir}/non_zero_region_contig.fasta
    rm ${out_dir}/non_zero_region.bed
    rm ${out_dir}/cord.txt 
    rm ${out_dir}/alignment.sam 
    rm ${out_dir}/alignment_sorted.bam
    rm ${out_dir}/alignment_sorted.bam.bai
    rm ${out_dir}/
    fi 

    

    if [[ -s $unpaired && $(wc -l < $unpaired) > 0 ]]; then

        # align reads
        #bowtie2 -x $index/index -U $unpaired --al ${out_dir}/reads_mapped.u.fq --un ${out_dir}/reads_unmapped.u.fq -S ${out_dir}/alignment.sam
        
        if [[ $check == 1 ]]; then
            # remove previously unmapped reads
           echo "testrun"
        else
          echo "testrun1"
        fi
        #unpaired=${out_dir}/reads_unmapped.u.fq
        # reassign unaligned reads
        

        tmp_num=$(wc -l < ${out_dir}/reads_mapped.u.fq)
        tmp_num=$((tmp_num/4 | bc))
        num_mapped_reads=$((num_mapped_reads + tmp_num))
        run=1
    fi

    # no reads left
    if [[ $run = 0 ]]; then
        break;
    fi
    run=0

done < $min_list

# log set of reads not mapped for de novo
if [[ -f $forward && -f $reverse ]]; then
    if [[ -f $unpaired ]]; then
        printf "%s\n" "--forward $forward --reverse $reverse --unpaired $unpaired" >> $log
    fi
    printf "%s\n" "--forward $forward --reverse $reverse" >> $log
fi

# check correct number of reads
if [[ -f $forward &&  -f $reverse ]]; then
    tmp_num=$(wc -l < $forward)
    num_unmapped_reads=$((tmp_num/4 | bc))
elif [[ -f $unpaired ]]; then
    tmp_num=$(wc -l < $unpaired)
    num_unmapped_reads=$((tmp_num/4 | bc))
fi
num_split_reads=$((num_mapped_reads + num_unmapped_reads))

if [[ $num_reads -ne $num_split_reads ]]; then

    printf "%s\n" "Metamorphic relation potentially violated: Number of original reads does not match number of mapped and unmapped reads!!!\n" >> $log
    printf "%s\n" "Please check parameters and contact support if unresolved!!!\n" >> $log
 
fi

