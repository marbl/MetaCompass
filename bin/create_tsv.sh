#!/bin/bash
#metacompass.final.ctg.fa 
final_contigs=$1
mc_contigs=$2
#mc_contigs=$2
#metacompass_output/metacompass.recruited.ids
recruited_ids=$3
output=$4

grep '>' $2 | tr -d '>'|sed 's/_[0-9]\+ / /g' >.tmp
cut -f1-2 -d '_' .tmp >.tmp_mc
grep '>' $final_contigs |tr -d '>' |grep ^k >.tmp_mg_ids
grep '>' $final_contigs |tr -d '>' |grep ^N >.tmp_mc_ids
samtools faidx $final_contigs |cut -f1-2>.tmp2
grep ^N "$final_contigs".fai|cut -f2 >.tmp2_mc_ln
grep ^k "$final_contigs".fai|cut -f2 >.tmp2_mg_ln

#if [ -f .tmp_tax ] ; then
#    rm .tmp_tax
#fi

if [ $recruited_ids !="NA" ] ; then
    for i in $(cat .tmp_mc);do
        grep $i $recruited_ids |cut -f2,4 >>.tmp_tax
    done
    printf "contig ID\tcontig size\treference genome\tposition start\tposition end\ttaxonomy id\tgenome name\n" >$output
    paste .tmp_mc_ids .tmp2_mc_ln .tmp .tmp_tax >>$output
else
    printf "contig ID\tcontig size\treference genome\tposition start\tposition end\n" >$output
    paste .tmp_mc_ids .tmp2_mc_ln .tmp  >>$output
fi

paste .tmp_mg_ids .tmp2_mg_ln >>$output
rm .tmp*

       
