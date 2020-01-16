#!/bin/bash
mc_contigs_pilon=$1 #contigs.pilon.fasta
mc_contigs=$2 #contigs.fasta
mg_contigs=$3 #megahit contigs
recruited_ids=$4 #metacompass_output/metacompass.recruited.ids or -r $file
assembled_fna=$5 #metacompass_output/metacompass.assembled.ids
output=$6 #output file
#echo $1 $2 $3 $4 $5
tmp=$RANDOM
#processing metacompass buildcontig output to get genome ids and genome start and end positions
if [ -s ${mc_contigs} ] ; then
	grep '>' ${mc_contigs} |rev| cut -f2- -d '_'|rev|tr -d '>'>.tmp${tmp}_genome_id	
	grep '>' ${mc_contigs} |rev| cut -f1,2 -d ' ' |rev >.tmp${tmp}_coords
fi

#processing pilon contigs to get ids and length
if [ -s ${mc_contigs_pilon} ] ; then
	grep '>' ${mc_contigs_pilon} |tr -d '>' >.tmp${tmp}_mc_ids	
	samtools faidx ${mc_contigs_pilon}
	cat ${mc_contigs_pilon}.fai|cut -f2 >.tmp${tmp}_mc_ln
fi

#processing megahit contigs to get ids and length
if [ -s ${mg_contigs} ] ; then
	grep '>' ${mg_contigs} |tr -d '>' >.tmp${tmp}_mg_ids
	samtools faidx ${mg_contigs}
	cat ${mg_contigs}.fai|cut -f2 >.tmp${tmp}_mg_ln
fi

if [ -s ${output} ] ; then
	rm $output
fi
#processing reference ids
if [ -s ${recruited_ids} ] ; then
    for i in $(cat .tmp${tmp}_genome_id);do
        grep $i $assembled_fna |cut -f2- -d ' ' >>.tmp${tmp}_name
    done
    printf "contig ID\tcontig size\treference genome\tposition start\tposition end\tgenome name\n" > $output
    paste .tmp${tmp}_mc_ids .tmp${tmp}_mc_ln .tmp${tmp}_genome_id .tmp${tmp}_coords .tmp${tmp}_name >> $output
#	rm .tmp${tmp}_genome_id .tmp${tmp}_coords .tmp${tmp}_name
elif [ -s ${mc_contigs_pilon} ] ; then
    	printf "contig ID\tcontig size\treference genome\tposition start\tposition end\n" > $output
    	paste .tmp${tmp}_mc_ids .tmp${tmp}_mc_ln >> $output
fi

#if [ -s .tmp${tmp}_mc_ids ];then
#	rm .tmp${tmp}_mc_ids
#fi
#if [ -s .tmp${tmp}_mc_ln ];then
#	rm .tmp${tmp}_mc_ln
#fi
 

#create final megahit +metacompass summary
if [ -s ${mg_contigs} ] ; then
	paste .tmp${tmp}_mg_ids .tmp${tmp}_mg_ln >> $output
#	rm .tmp${tmp}_mg_ids .tmp${tmp}_mg_ln
fi
rm .tmp${tmp}*

       
